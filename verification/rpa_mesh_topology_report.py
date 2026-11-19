#!/usr/bin/env python3
"""Geometry-only TRI/QUAD audit. h_axi is NOT a measured physical CFL step."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import xml.etree.ElementTree as ET


def percentile(values, fraction):
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    left = int(position)
    return ordered[left] + (ordered[min(left + 1, len(ordered)-1)]-ordered[left])*(position-left)


def audit(path):
    piece = ET.parse(path).find('.//Piece')
    flat = list(map(float, piece.find('Points/DataArray').text.split()))
    points = list(zip(flat[::3], flat[1::3]))
    data = {a.attrib['Name']: list(map(float, a.text.split())) for a in piece.findall('CellData/DataArray')}
    topology = {a.attrib['Name']: list(map(int, a.text.split())) for a in piece.findall('Cells/DataArray')}
    cells, faces, start = [], {}, 0
    for index, (end, kind) in enumerate(zip(topology['offsets'], topology['types'])):
        ids = topology['connectivity'][start:end]
        assert kind in (5, 9) and len(ids) == (3 if kind == 5 else 4)
        polygon = [points[n] for n in ids]
        edges = list(zip(polygon, polygon[1:]+polygon[:1]))
        crosses = [a[0]*b[1]-b[0]*a[1] for a, b in edges]
        twice_area = sum(crosses)
        assert twice_area > 0
        area = .5*twice_area
        centre = [sum((a[k]+b[k])*cross for (a,b),cross in zip(edges,crosses))/(3*twice_area) for k in (0,1)]
        lengths = [math.dist(a,b) for a,b in edges]
        surface = sum(2*math.pi*.5*(a[1]+b[1])*length for (a,b),length in zip(edges,lengths))
        volume = data['axisymmetric_volume'][index]
        assert all(math.isfinite(value) and value > 0 for value in (area,volume,surface))
        angles=[]
        for i in range(len(polygon)):
            a,b,c=polygon[i],polygon[(i+1)%len(polygon)],polygon[(i+2)%len(polygon)]
            angles.append(((b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]))/(math.dist(a,b)*math.dist(b,c)))
        row = dict(cell=index,region=int(data['rpa_region'][index]),type=kind,x=centre[0],r=centre[1],
                   area=area,volume=volume,edge_min=min(lengths),edge_max=max(lengths),
                   dx=max(p[0] for p in polygon)-min(p[0] for p in polygon),
                   dr=max(p[1] for p in polygon)-min(p[1] for p in polygon),
                   size=math.sqrt(area),quality=data['cell_quality'][index],
                   aspect=max(lengths)**2/(area if kind==9 else 2*area),h_axi=volume/surface,
                   angular_quality=min(angles),block=int(data.get('rpa_block',data['rpa_region'])[index]),
                   nonorthogonality=0.0,gradation=1.0,polygon=polygon)
        cells.append(row)
        for a,b in zip(ids,ids[1:]+ids[:1]):
            faces.setdefault(tuple(sorted((a,b))),[]).append(index)
        start=end
    for (a,b), owners in faces.items():
        assert len(owners) in (1,2), 'Nonmanifold face'
        p,q=points[a],points[b]
        normal=(q[1]-p[1],p[0]-q[0])
        c=cells[owners[0]]
        target=([cells[owners[1]]['x'],cells[owners[1]]['r']] if len(owners)==2 else [(p[0]+q[0])/2,(p[1]+q[1])/2])
        d=(target[0]-c['x'],target[1]-c['r'])
        cosine=abs(normal[0]*d[0]+normal[1]*d[1])/(math.hypot(*normal)*math.hypot(*d))
        angle=math.degrees(math.acos(min(1.0,cosine)))
        for owner in owners:
            cells[owner]['nonorthogonality']=max(cells[owner]['nonorthogonality'],angle)
            if len(owners)==2:
                ratio=cells[owners[0]]['size']/cells[owners[1]]['size']
                cells[owner]['gradation']=max(cells[owner]['gradation'],ratio,1/ratio)
    regions={}
    for region in range(1,6):
        group=[c for c in cells if c['region']==region]
        assert group
        stats=dict(cells=len(group),tri=sum(c['type']==5 for c in group),quad=sum(c['type']==9 for c in group))
        for name in ('edge_min','edge_max','dx','dr','area','volume','quality','angular_quality','aspect','h_axi','size','nonorthogonality','gradation'):
            values=[c[name] for c in group]
            stats[name]=dict(min=min(values),max=max(values),mean=statistics.fmean(values),
                             p01=percentile(values,.01),p50=percentile(values,.5),p99=percentile(values,.99))
        regions[region]=stats
    return dict(cells=len(cells),faces=len(faces),regions=regions,
                smallest_h_axi=[{k:v for k,v in c.items() if k!='polygon'} for c in sorted(cells,key=lambda c:c['h_axi'])[:20]]),cells


def svg(path, cells, limits):
    left,right,bottom,top=limits
    width=1200; height=max(240, min(800, width*(top-bottom)/(right-left)))
    palette=['#dddddd','#749fda','#e5a76f','#94c997','#c8df83','#69c9bb','#b2c5c8','#99b0c8','#d7b4df']
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="white"/>']
    for cell in cells:
        polygon=cell['polygon']
        if max(p[0] for p in polygon)<left or min(p[0] for p in polygon)>right or max(p[1] for p in polygon)<bottom or min(p[1] for p in polygon)>top:continue
        coords=' '.join(f'{(x-left)/(right-left)*width:.5f},{height-(y-bottom)/(top-bottom)*height:.5f}' for x,y in polygon)
        lines.append(f'<polygon points="{coords}" fill="{palette[cell["block"]]}" stroke="#333" stroke-width=".25"/>')
    lines.append('</svg>')
    path.write_text('\n'.join(lines))


def write_subsets(source, cells, output):
    """Visualization subsets only; retain original field tokens/precision."""
    piece=ET.parse(source).find('.//Piece')
    xyz=piece.find('Points/DataArray').text.split()
    topology={a.attrib['Name']:list(map(int,a.text.split())) for a in piece.findall('Cells/DataArray')}
    offsets=[0]+topology['offsets']
    for block,name in ((1,'engine'),(3,'plume_block_1'),(4,'plume_block_2'),(5,'plume_block_3'),(8,'transition')):
        selected=[c['cell'] for c in cells if c['block']==block]
        if not selected:continue
        source_nodes=sorted({n for i in selected for n in topology['connectivity'][offsets[i]:offsets[i+1]]})
        mapping={node:i for i,node in enumerate(source_nodes)}
        root=ET.Element('VTKFile',type='UnstructuredGrid',version='0.1',byte_order='LittleEndian')
        target=ET.SubElement(ET.SubElement(root,'UnstructuredGrid'),'Piece',NumberOfPoints=str(len(source_nodes)),NumberOfCells=str(len(selected)))
        points=ET.SubElement(target,'Points')
        ET.SubElement(points,'DataArray',type='Float64',NumberOfComponents='3',format='ascii').text=' '.join(v for node in source_nodes for v in xyz[3*node:3*node+3])
        connectivity=[]; ends=[]
        for i in selected:
            connectivity.extend(mapping[n] for n in topology['connectivity'][offsets[i]:offsets[i+1]])
            ends.append(len(connectivity))
        dest=ET.SubElement(target,'Cells')
        for key,kind,values in (('connectivity','Int64',connectivity),('offsets','Int64',ends),('types','UInt8',[topology['types'][i] for i in selected])):
            ET.SubElement(dest,'DataArray',Name=key,type=kind,format='ascii').text=' '.join(map(str,values))
        dest=ET.SubElement(target,'CellData')
        for array in piece.findall('CellData/DataArray'):
            count=int(array.attrib.get('NumberOfComponents',1));values=array.text.split()
            ET.SubElement(dest,'DataArray',array.attrib).text=' '.join(v for i in selected for v in values[i*count:(i+1)*count])
        ET.ElementTree(root).write(output/(name+'.vtu'),encoding='utf-8',xml_declaration=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mesh',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--subsets',action='store_true',help='Export standalone ENGINE/PLUME-block/TRANSITION VTUs')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    report,cells=audit(args.mesh)
    (args.output/'metrics.json').write_text(json.dumps(report,indent=2)+'\n')
    with (args.output/'cells.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=[key for key in cells[0] if key!='polygon'])
        writer.writeheader();writer.writerows({k:v for k,v in cell.items() if k!='polygon'} for cell in cells)
    xmax=max(p[0] for cell in cells for p in cell['polygon']);rmax=max(p[1] for cell in cells for p in cell['polygon'])
    svg(args.output/'global.svg',cells,(0,xmax,0,rmax))
    engine=[c for c in cells if c['region']==1];exit_x=max(p[0] for c in engine for p in c['polygon']);exit_r=max(c['r'] for c in engine if c['x']>.99*exit_x)
    svg(args.output/'engine-and-exit.svg',cells,(0,exit_x+4*exit_r,0,4*exit_r))
    if args.subsets:write_subsets(args.mesh,cells,args.output)
    print(args.output,report['cells'],'cells',report['faces'],'faces; h_axi is a geometric proxy only')


if __name__=='__main__':main()
