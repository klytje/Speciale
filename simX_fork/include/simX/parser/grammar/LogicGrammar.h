//
// Created by munk on 13-12-15.
//

#ifndef SIMX_PARSER_LOGICGRAMMAR_H
#define SIMX_PARSER_LOGICGRAMMAR_H

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace simX {
    namespace parser {
        namespace logic {
            namespace qi = boost::spirit::qi;
            namespace ascii = boost::spirit::ascii;

            typedef std::string::const_iterator Iterator;


            template <typename tag, typename T> struct binop;
            template <typename tag, typename T> struct unop;

            struct op_or  {};
            struct op_and {};
            struct op_xor {};
            struct op_not {};

            template <typename T>
            using LogicExpression = boost::variant<T,
                    boost::recursive_wrapper<unop <op_not, T> >,
                    boost::recursive_wrapper<binop<op_and, T> >,
                    boost::recursive_wrapper<binop<op_xor, T> >,
                    boost::recursive_wrapper<binop<op_or, T> >
            >;

            template <typename tag, typename T> struct binop
            {
                explicit binop(const LogicExpression<T>& l, const LogicExpression<T>& r) : oper1(l), oper2(r) { }
                LogicExpression<T> oper1, oper2;
            };

            template <typename tag, typename T> struct unop
            {
                explicit unop(const LogicExpression<T>& o) : oper1(o) { }
                LogicExpression<T> oper1;
            };



            template <typename T, typename It, typename Skipper = qi::space_type>
            struct LogicGrammar : qi::grammar<It, LogicExpression<T>(), Skipper>
            {
                LogicGrammar() : LogicGrammar::base_type(expr_)
                {
                    using namespace qi;
                    namespace phx   = boost::phoenix;
                    using boost::spirit::ascii::alnum;

                    expr_  = or_.alias();

                    or_  = (xor_ >> "|"  >> or_ ) [ qi::_val = phx::construct<binop<op_or,T>>(qi::_1, qi::_2) ] | xor_   [ qi::_val = qi::_1 ];
                    xor_ = (and_ >> "^" >> xor_) [ qi::_val = phx::construct<binop<op_xor,T>>(qi::_1, qi::_2) ] | and_   [ qi::_val = qi::_1 ];
                    and_ = (not_ >> "&" >> and_) [ qi::_val = phx::construct<binop<op_and,T>>(qi::_1, qi::_2) ] | not_   [ qi::_val = qi::_1 ];
                    not_ = ("!" > simple       ) [ qi::_val = phx::construct<unop <op_not,T>>(qi::_1)     ] | simple [ qi::_val = qi::_1 ];

                    simple = (('(' > expr_ > ')') | var_);
                    var_ = qi::lexeme[ +alnum ] >> eps;
                }

            private:
                qi::rule<It, T() , Skipper> var_;
                qi::rule<It, LogicExpression<T>(), Skipper> not_, and_, xor_, or_, simple, expr_;
            };

            template <typename T, typename It, typename Skipper = qi::space_type>
            struct DownscaledGrammar : qi::grammar<It, LogicExpression<T>(), Skipper>
            {
                DownscaledGrammar() : DownscaledGrammar::base_type(expr_)
                {
                    using namespace qi;
                    namespace phx   = boost::phoenix;
                    using boost::spirit::ascii::alnum;
                    using boost::spirit::uint_;
                    using qi::lit;

                    expr_  = or_.alias();

                    or_  = (xor_ >> "|" >> -lit('|') > or_ ) [ qi::_val = phx::construct<binop<op_or,T>>(qi::_1, qi::_2) ] | xor_   [ qi::_val = qi::_1 ];
                    xor_ = (and_ >> "^" > xor_) [ qi::_val = phx::construct<binop<op_xor,T>>(qi::_1, qi::_2) ] | and_   [ qi::_val = qi::_1 ];
                    and_ = (not_ >> "&" >> -lit('&')  > and_) [ qi::_val = phx::construct<binop<op_and,T>>(qi::_1, qi::_2) ] | not_   [ qi::_val = qi::_1 ];
                    not_ = ("!" > simple       ) [ qi::_val = phx::construct<unop <op_not,T>>(qi::_1)     ] | simple [ qi::_val = qi::_1 ];

                    simple = (('(' > expr_ > ')') | var_);
                    var_ = qi::lexeme[ +alnum ] >> (('(' > uint_ > ')') | qi::attr(1U)) >> eps;

                    or_.name("OR");
                    and_.name("AND");
                    xor_.name("XOR");
                    expr_.name("expression");
                    not_.name("expression");
                    simple.name("simple");
                    var_.name("variable");

                    on_error<fail>(expr_,
                                   phx::ref(std::cerr)
                                           << "Error! Expecting "
                                           << qi::_4
                                           << " got: '"
                                           << phx::construct<std::string>(qi::_3, qi::_2)
                                           << "'\n"
                    );
                }

            private:
                qi::rule<It, T() , Skipper> var_;
                qi::rule<It, LogicExpression<T>(), Skipper> not_, and_, xor_, or_, simple, expr_;
            };
        }
    }
}
#endif //SIMX_PARSER_LOGICGRAMMAR_H
