#include "cfd/input/OpenFOAMCaseReader.hpp"

#include "cfd/mesh/Boundary.hpp"
#include "cfd/mesh/Cell.hpp"
#include "cfd/mesh/Mesh.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace cfd::input
{
namespace
{

enum class TokenKind : std::uint8_t
{
    Word,
    Number,
    QuotedString,
    LeftBrace,
    RightBrace,
    LeftBracket,
    RightBracket,
    Semicolon,
    End
};

struct Token
{
    TokenKind kind{TokenKind::End};
    std::string_view text;
    std::size_t line{1};
    double number{};
};

[[nodiscard]]
std::string read_file(const std::filesystem::path &file_path)
{
    std::ifstream input{file_path, std::ios::binary};

    if (!input)
    {
        throw std::runtime_error(file_path.string() + ": Unable to open input file.");
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();

    if (!input.eof() && input.fail())
    {
        throw std::runtime_error(file_path.string() + ": Error while reading input file.");
    }

    return buffer.str();
}

class CaseParser
{
  public:
    explicit CaseParser(std::filesystem::path file_path)
        : file_path_(std::move(file_path)), source_(read_file(file_path_))
    {
        advance();
    }

    [[nodiscard]]
    MeshInput parse_mesh_dict()
    {
        static_cast<void>(parse_foam_header("dictionary", "meshDict"));

        const RectangleGeometry geometry{parse_geometry()};
        const MeshGenerationOptions generation_options{parse_mesh_options()};

        expect(TokenKind::End, "end of file");

        return {
            .geometry = geometry,
            .generation_options = generation_options,
        };
    }

    [[nodiscard]]
    ScalarFieldInput parse_scalar_field()
    {
        std::string object_name{parse_foam_header("volScalarField", {})};

        expect_keyword("dimensions");
        expect(TokenKind::LeftBracket, "'['");

        constexpr std::size_t dimension_count{7};
        for (std::size_t dimension_index = 0; dimension_index < dimension_count; ++dimension_index)
        {
            const Token dimension{take_number("a dimension exponent")};

            if (dimension.number != 0.0)
            {
                fail(dimension.line, "Only dimensionless scalar fields are supported.");
            }
        }

        expect(TokenKind::RightBracket, "']'");
        expect(TokenKind::Semicolon, "';'");

        expect_keyword("internalField");
        const double internal_value{parse_uniform_scalar("internalField")};
        expect(TokenKind::Semicolon, "';'");

        expect_keyword("boundaryField");
        expect(TokenKind::LeftBrace, "'{'");

        std::vector<NamedScalarBoundaryCondition> boundary_conditions;

        while (current_.kind != TokenKind::RightBrace)
        {
            if (current_.kind == TokenKind::End)
            {
                fail_current("Expected a boundary entry or '}'.");
            }

            const Token boundary_name_token{take_word("a boundary name")};
            const std::string boundary_name{boundary_name_token.text};

            const auto duplicate{std::ranges::find_if(
                boundary_conditions, [&boundary_name](const NamedScalarBoundaryCondition &named_condition) {
                    return named_condition.boundary_name == boundary_name;
                })};

            if (duplicate != boundary_conditions.end())
            {
                fail(boundary_name_token.line, "Duplicate boundary entry '" + boundary_name + "'.");
            }

            expect(TokenKind::LeftBrace, "'{'");
            ScalarBoundaryCondition condition{parse_boundary_condition(boundary_name)};
            expect(TokenKind::RightBrace, "'}'");

            boundary_conditions.push_back({
                .boundary_name = boundary_name,
                .condition = condition,
            });
        }

        expect(TokenKind::RightBrace, "'}'");
        expect(TokenKind::End, "end of file");

        return {
            .object_name = std::move(object_name),
            .internal_value = internal_value,
            .boundary_conditions = std::move(boundary_conditions),
        };
    }

  private:
    [[noreturn]]
    void fail(const std::size_t line, const std::string &message) const
    {
        throw std::runtime_error(file_path_.string() + ':' + std::to_string(line) + ": " + message);
    }

    [[noreturn]]
    void fail_current(const std::string &message) const
    {
        fail(current_.line, message);
    }

    [[nodiscard]]
    static bool is_word_start(const char character) noexcept
    {
        const auto value{static_cast<unsigned char>(character)};
        return std::isalpha(value) != 0 || character == '_';
    }

    [[nodiscard]]
    static bool is_word_character(const char character) noexcept
    {
        const auto value{static_cast<unsigned char>(character)};
        return std::isalnum(value) != 0 || character == '_' || character == '-' || character == '.';
    }

    [[nodiscard]]
    bool is_number_start() const noexcept
    {
        const char character{source_[position_]};

        if (std::isdigit(static_cast<unsigned char>(character)) != 0)
        {
            return true;
        }

        if (character == '.')
        {
            return position_ + 1 < source_.size() &&
                   std::isdigit(static_cast<unsigned char>(source_[position_ + 1])) != 0;
        }

        if (character != '+' && character != '-')
        {
            return false;
        }

        if (position_ + 1 >= source_.size())
        {
            return false;
        }

        const char next{source_[position_ + 1]};
        return std::isdigit(static_cast<unsigned char>(next)) != 0 ||
               (next == '.' && position_ + 2 < source_.size() &&
                std::isdigit(static_cast<unsigned char>(source_[position_ + 2])) != 0);
    }

    [[nodiscard]]
    bool at_end() const noexcept
    {
        return position_ >= source_.size();
    }

    void skip_ignored()
    {
        while (!at_end())
        {
            const char character{source_[position_]};

            if (std::isspace(static_cast<unsigned char>(character)) != 0)
            {
                if (character == '\n')
                {
                    ++line_;
                }

                ++position_;
                continue;
            }

            if (character != '/' || position_ + 1 >= source_.size())
            {
                return;
            }

            const char next{source_[position_ + 1]};

            if (next == '/')
            {
                position_ += 2;

                while (!at_end() && source_[position_] != '\n')
                {
                    ++position_;
                }

                continue;
            }

            if (next != '*')
            {
                return;
            }

            const std::size_t comment_line{line_};
            position_ += 2;

            while (position_ + 1 < source_.size() && !(source_[position_] == '*' && source_[position_ + 1] == '/'))
            {
                if (source_[position_] == '\n')
                {
                    ++line_;
                }

                ++position_;
            }

            if (position_ + 1 >= source_.size())
            {
                fail(comment_line, "Unterminated block comment.");
            }

            position_ += 2;
        }
    }

    [[nodiscard]]
    Token lex_number()
    {
        const std::size_t start{position_};
        const std::size_t token_line{line_};

        if (source_[position_] == '+' || source_[position_] == '-')
        {
            ++position_;
        }

        while (!at_end() && std::isdigit(static_cast<unsigned char>(source_[position_])) != 0)
        {
            ++position_;
        }

        if (!at_end() && source_[position_] == '.')
        {
            ++position_;

            while (!at_end() && std::isdigit(static_cast<unsigned char>(source_[position_])) != 0)
            {
                ++position_;
            }
        }

        if (!at_end() && (source_[position_] == 'e' || source_[position_] == 'E'))
        {
            ++position_;

            if (!at_end() && (source_[position_] == '+' || source_[position_] == '-'))
            {
                ++position_;
            }

            if (at_end() || std::isdigit(static_cast<unsigned char>(source_[position_])) == 0)
            {
                fail(token_line, "Malformed scalar exponent.");
            }

            while (!at_end() && std::isdigit(static_cast<unsigned char>(source_[position_])) != 0)
            {
                ++position_;
            }
        }

        if (!at_end() && is_word_character(source_[position_]))
        {
            fail(token_line, "Malformed scalar token.");
        }

        const std::string_view text{std::string_view{source_}.substr(start, position_ - start)};
        const std::string text_copy{text};
        std::size_t converted_character_count{};
        double value{};

        try
        {
            value = std::stod(text_copy, &converted_character_count);
        }
        catch (const std::invalid_argument &)
        {
            fail(token_line, "Scalar values must be finite decimal numbers.");
        }
        catch (const std::out_of_range &)
        {
            fail(token_line, "Scalar values must be finite decimal numbers.");
        }

        if (converted_character_count != text_copy.size() || !std::isfinite(value))
        {
            fail(token_line, "Scalar values must be finite decimal numbers.");
        }

        return {
            .kind = TokenKind::Number,
            .text = text,
            .line = token_line,
            .number = value,
        };
    }

    [[nodiscard]]
    Token lex_quoted_string()
    {
        const std::size_t token_line{line_};
        ++position_;
        const std::size_t start{position_};

        while (!at_end() && source_[position_] != '"')
        {
            if (source_[position_] == '\n')
            {
                fail(token_line, "Quoted strings may not span lines.");
            }

            if (source_[position_] == '\\')
            {
                fail(token_line, "Escape sequences in quoted strings are unsupported.");
            }

            ++position_;
        }

        if (at_end())
        {
            fail(token_line, "Unterminated quoted string.");
        }

        const std::string_view text{std::string_view{source_}.substr(start, position_ - start)};
        ++position_;

        return {
            .kind = TokenKind::QuotedString,
            .text = text,
            .line = token_line,
        };
    }

    [[nodiscard]]
    Token lex_token()
    {
        skip_ignored();

        if (at_end())
        {
            return {
                .kind = TokenKind::End,
                .text = {},
                .line = line_,
            };
        }

        const std::size_t token_line{line_};
        const char character{source_[position_]};

        switch (character)
        {
        case '{':
            ++position_;
            return {.kind = TokenKind::LeftBrace, .text = "{", .line = token_line};
        case '}':
            ++position_;
            return {.kind = TokenKind::RightBrace, .text = "}", .line = token_line};
        case '[':
            ++position_;
            return {.kind = TokenKind::LeftBracket, .text = "[", .line = token_line};
        case ']':
            ++position_;
            return {.kind = TokenKind::RightBracket, .text = "]", .line = token_line};
        case ';':
            ++position_;
            return {.kind = TokenKind::Semicolon, .text = ";", .line = token_line};
        case '"':
            return lex_quoted_string();
        case '#':
            fail(token_line, "Preprocessing directives are unsupported.");
        case '$':
            fail(token_line, "Variable substitution is unsupported.");
        case '(':
        case ')':
            fail(token_line, "Vector values are unsupported by the scalar-field reader.");
        default:
            break;
        }

        if (is_number_start())
        {
            return lex_number();
        }

        if (!is_word_start(character))
        {
            fail(token_line, "Unsupported character '" + std::string(1, character) + "'.");
        }

        const std::size_t start{position_};
        ++position_;

        while (!at_end() && is_word_character(source_[position_]))
        {
            ++position_;
        }

        return {
            .kind = TokenKind::Word,
            .text = std::string_view{source_}.substr(start, position_ - start),
            .line = token_line,
        };
    }

    void advance()
    {
        current_ = lex_token();
    }

    void expect(const TokenKind kind, const std::string_view description)
    {
        if (current_.kind != kind)
        {
            fail_current("Expected " + std::string(description) + ".");
        }

        advance();
    }

    void expect_keyword(const std::string_view keyword)
    {
        if (current_.kind != TokenKind::Word || current_.text != keyword)
        {
            fail_current("Expected '" + std::string(keyword) + "'.");
        }

        advance();
    }

    [[nodiscard]]
    Token take_word(const std::string_view description)
    {
        if (current_.kind != TokenKind::Word && current_.kind != TokenKind::QuotedString)
        {
            fail_current("Expected " + std::string(description) + ".");
        }

        const Token result{current_};
        advance();
        return result;
    }

    [[nodiscard]]
    Token take_key(const std::string_view dictionary_name)
    {
        if (current_.kind != TokenKind::Word)
        {
            fail_current("Expected an entry name in " + std::string(dictionary_name) + ".");
        }

        const Token result{current_};
        advance();
        return result;
    }

    [[nodiscard]]
    Token take_number(const std::string_view description)
    {
        if (current_.kind != TokenKind::Number)
        {
            fail_current("Expected " + std::string(description) + ".");
        }

        const Token result{current_};
        advance();
        return result;
    }

    [[nodiscard]]
    double parse_uniform_scalar(const std::string_view entry_name)
    {
        const Token distribution{take_word("'uniform'")};

        if (distribution.text == "nonuniform")
        {
            fail(distribution.line, "Nonuniform " + std::string(entry_name) + " values are unsupported.");
        }

        if (distribution.text != "uniform")
        {
            fail(distribution.line, "Expected 'uniform' for " + std::string(entry_name) + ".");
        }

        return take_number("a finite scalar value").number;
    }

    [[nodiscard]]
    std::string parse_foam_header(const std::string_view expected_class, const std::string_view expected_object)
    {
        expect_keyword("FoamFile");
        expect(TokenKind::LeftBrace, "'{'");

        std::optional<Token> version;
        std::optional<Token> format;
        std::optional<Token> class_name;
        std::optional<Token> object_name;
        std::optional<Token> location;

        while (current_.kind != TokenKind::RightBrace)
        {
            const Token key{take_key("FoamFile")};

            if (key.text == "version")
            {
                if (version.has_value())
                {
                    fail(key.line, "Duplicate FoamFile version entry.");
                }

                version = take_number("the FoamFile version");
            }
            else if (key.text == "format")
            {
                if (format.has_value())
                {
                    fail(key.line, "Duplicate FoamFile format entry.");
                }

                format = take_word("the FoamFile format");
            }
            else if (key.text == "class")
            {
                if (class_name.has_value())
                {
                    fail(key.line, "Duplicate FoamFile class entry.");
                }

                class_name = take_word("the FoamFile class");
            }
            else if (key.text == "object")
            {
                if (object_name.has_value())
                {
                    fail(key.line, "Duplicate FoamFile object entry.");
                }

                object_name = take_word("the FoamFile object");
            }
            else if (key.text == "location")
            {
                if (location.has_value())
                {
                    fail(key.line, "Duplicate FoamFile location entry.");
                }

                location = take_word("the FoamFile location");
            }
            else
            {
                fail(key.line, "Unsupported FoamFile entry '" + std::string(key.text) + "'.");
            }

            expect(TokenKind::Semicolon, "';'");
        }

        expect(TokenKind::RightBrace, "'}'");

        if (!format.has_value() || !class_name.has_value() || !object_name.has_value())
        {
            fail_current("FoamFile requires format, class, and object entries.");
        }

        if (version.has_value() && version->number != 2.0)
        {
            fail(version->line, "Only FoamFile version 2.0 is supported.");
        }

        if (format->text != "ascii")
        {
            fail(format->line, "Only FoamFile format 'ascii' is supported.");
        }

        if (class_name->text != expected_class)
        {
            fail(class_name->line, "FoamFile class must be '" + std::string(expected_class) + "'.");
        }

        if (!expected_object.empty() && object_name->text != expected_object)
        {
            fail(object_name->line, "FoamFile object must be '" + std::string(expected_object) + "'.");
        }

        return std::string{object_name->text};
    }

    [[nodiscard]]
    RectangleGeometry parse_geometry()
    {
        expect_keyword("geometry");
        expect(TokenKind::LeftBrace, "'{'");

        std::optional<Token> type;
        std::optional<Token> length;
        std::optional<Token> height;

        while (current_.kind != TokenKind::RightBrace)
        {
            const Token key{take_key("geometry")};

            if (key.text == "type")
            {
                if (type.has_value())
                {
                    fail(key.line, "Duplicate geometry type entry.");
                }

                type = take_word("the geometry type");
            }
            else if (key.text == "length")
            {
                if (length.has_value())
                {
                    fail(key.line, "Duplicate geometry length entry.");
                }

                length = take_number("a finite rectangle length");
            }
            else if (key.text == "height")
            {
                if (height.has_value())
                {
                    fail(key.line, "Duplicate geometry height entry.");
                }

                height = take_number("a finite rectangle height");
            }
            else
            {
                fail(key.line, "Unsupported geometry entry '" + std::string(key.text) + "'.");
            }

            expect(TokenKind::Semicolon, "';'");
        }

        expect(TokenKind::RightBrace, "'}'");

        if (!type.has_value() || !length.has_value() || !height.has_value())
        {
            fail_current("geometry requires type, length, and height entries.");
        }

        if (type->text != "rectangle")
        {
            fail(type->line, "Only rectangle geometry is supported.");
        }

        if (length->number <= 0.0)
        {
            fail(length->line, "Rectangle length must be finite and positive.");
        }

        if (height->number <= 0.0)
        {
            fail(height->line, "Rectangle height must be finite and positive.");
        }

        return {
            .length = length->number,
            .height = height->number,
        };
    }

    [[nodiscard]]
    MeshGenerationOptions parse_mesh_options()
    {
        expect_keyword("mesh");
        expect(TokenKind::LeftBrace, "'{'");

        std::optional<Token> cell_type_name;
        std::optional<Token> mesh_size;

        while (current_.kind != TokenKind::RightBrace)
        {
            const Token key{take_key("mesh")};

            if (key.text == "cellType")
            {
                if (cell_type_name.has_value())
                {
                    fail(key.line, "Duplicate mesh cellType entry.");
                }

                cell_type_name = take_word("the mesh cellType");
            }
            else if (key.text == "size")
            {
                if (mesh_size.has_value())
                {
                    fail(key.line, "Duplicate mesh size entry.");
                }

                mesh_size = take_number("a finite mesh size");
            }
            else
            {
                fail(key.line, "Unsupported mesh entry '" + std::string(key.text) + "'.");
            }

            expect(TokenKind::Semicolon, "';'");
        }

        expect(TokenKind::RightBrace, "'}'");

        if (!cell_type_name.has_value() || !mesh_size.has_value())
        {
            fail_current("mesh requires cellType and size entries.");
        }

        CellType cell_type{};
        if (cell_type_name->text == "triangle")
        {
            cell_type = CellType::Triangle;
        }
        else if (cell_type_name->text == "quadrilateral")
        {
            cell_type = CellType::Quadrilateral;
        }
        else
        {
            fail(cell_type_name->line, "Unsupported mesh cellType '" + std::string(cell_type_name->text) + "'.");
        }

        if (mesh_size->number <= 0.0)
        {
            fail(mesh_size->line, "Mesh size must be finite and positive.");
        }

        return {
            .mesh_size = mesh_size->number,
            .cell_type = cell_type,
        };
    }

    [[nodiscard]]
    ScalarBoundaryCondition parse_boundary_condition(const std::string_view boundary_name)
    {
        std::optional<Token> type;
        std::optional<double> value;
        std::optional<double> gradient;

        while (current_.kind != TokenKind::RightBrace)
        {
            const Token key{take_key("boundary '" + std::string(boundary_name) + "'")};

            if (key.text == "type")
            {
                if (type.has_value())
                {
                    fail(key.line, "Duplicate boundary-condition type entry.");
                }

                type = take_word("the boundary-condition type");
            }
            else if (key.text == "value")
            {
                if (value.has_value())
                {
                    fail(key.line, "Duplicate boundary-condition value entry.");
                }

                value = parse_uniform_scalar("boundary value");
            }
            else if (key.text == "gradient")
            {
                if (gradient.has_value())
                {
                    fail(key.line, "Duplicate boundary-condition gradient entry.");
                }

                gradient = parse_uniform_scalar("boundary gradient");
            }
            else
            {
                fail(key.line, "Unsupported boundary-condition entry '" + std::string(key.text) + "'.");
            }

            expect(TokenKind::Semicolon, "';'");
        }

        if (!type.has_value())
        {
            fail_current("Boundary '" + std::string(boundary_name) + "' requires a type entry.");
        }

        if (type->text == "fixedValue")
        {
            if (!value.has_value() || gradient.has_value())
            {
                fail_current("fixedValue requires exactly one uniform value entry.");
            }

            return {ScalarBoundaryConditionType::Dirichlet, *value};
        }

        if (type->text == "zeroGradient")
        {
            if (value.has_value() || gradient.has_value())
            {
                fail_current("zeroGradient does not accept value or gradient entries.");
            }

            return {ScalarBoundaryConditionType::Neumann, 0.0};
        }

        if (type->text == "fixedGradient")
        {
            if (value.has_value() || !gradient.has_value())
            {
                fail_current("fixedGradient requires exactly one uniform gradient entry.");
            }

            return {ScalarBoundaryConditionType::Neumann, *gradient};
        }

        fail(type->line, "Unsupported scalar boundary-condition type '" + std::string(type->text) + "' on boundary '" +
                             std::string(boundary_name) + "'.");
    }

    std::filesystem::path file_path_;
    std::string source_;
    std::size_t position_{};
    std::size_t line_{1};
    Token current_;
};

} // namespace

MeshInput read_mesh_dict(const std::filesystem::path &file_path)
{
    CaseParser parser{file_path};
    return parser.parse_mesh_dict();
}

ScalarFieldInput read_scalar_field(const std::filesystem::path &file_path)
{
    CaseParser parser{file_path};
    return parser.parse_scalar_field();
}

ScalarBoundaryConditions resolve_boundary_conditions(const Mesh &mesh, const ScalarFieldInput &field_input)
{
    const std::span<const BoundaryGroup> boundary_groups{mesh.boundary_groups()};
    std::vector<const ScalarBoundaryCondition *> conditions_by_id(boundary_groups.size(), nullptr);

    for (const NamedScalarBoundaryCondition &named_condition : field_input.boundary_conditions)
    {
        const auto group{std::ranges::find_if(boundary_groups, [&named_condition](const BoundaryGroup &candidate) {
            return candidate.name == named_condition.boundary_name;
        })};

        if (group == boundary_groups.end())
        {
            throw std::invalid_argument("Unknown Mesh boundary '" + named_condition.boundary_name + "'.");
        }

        if (conditions_by_id[group->id] != nullptr)
        {
            throw std::invalid_argument("Duplicate scalar condition for boundary '" + named_condition.boundary_name +
                                        "'.");
        }

        conditions_by_id[group->id] = &named_condition.condition;
    }

    std::vector<ScalarBoundaryCondition> ordered_conditions;
    ordered_conditions.reserve(boundary_groups.size());

    for (const BoundaryGroup &group : boundary_groups)
    {
        const ScalarBoundaryCondition *const condition{conditions_by_id[group.id]};

        if (condition == nullptr)
        {
            throw std::invalid_argument("Missing scalar condition for boundary '" + group.name + "'.");
        }

        ordered_conditions.push_back(*condition);
    }

    return ScalarBoundaryConditions{boundary_groups.size(), std::move(ordered_conditions)};
}

} // namespace cfd::input
