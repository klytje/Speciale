//
// Created by munk on 13-12-15.
//

#include "simX/parser/DAQTriggerParser.h"
#include "simX/parser/grammar/LogicGrammar.h"
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#include <spdlog/details/format.h>
#include <simX/Logger.h>

using namespace std;
using namespace simX::parser::logic;
using simX::daq::SimpleDAQ;
using simX::detection::DetectionSystem;

namespace {
    struct Indexed {
        std::string name;
        size_t reduction;
        size_t count;
        size_t index;
    };

    struct InitStructures : boost::static_visitor<void>
    {
        InitStructures(const vector<string>& s) : names(s) {}
        const vector<string>& names;

        typedef Indexed var;
        typedef LogicExpression<Indexed> expr;

        //
        void operator()(var& v) const {
            v.count = 0;
            for (size_t i = 0; i < names.size(); ++i) {
                if (v.name == names[i]) {
                    v.index = i;
                    return;
                }
            }
            throw std::invalid_argument("There is no detector named '" + v.name + "'");
        }

        void operator()(binop<op_and, Indexed>& b) const { print(b.oper1, b.oper2); }
        void operator()(binop<op_or, Indexed>& b) const { print(b.oper1, b.oper2); }
        void operator()(binop<op_xor, Indexed>& b) const { print(b.oper1, b.oper2); }

        void print(expr& l, expr& r) const {
            boost::apply_visitor(*this, l);
            boost::apply_visitor(*this, r);
        }

        void operator()(unop<op_not, Indexed>& u) const {
            boost::apply_visitor(*this, u.oper1);
        }
    };

    struct IncrementCounts : boost::static_visitor<void>
    {
        IncrementCounts(const SimpleDAQ::TriggerStatus& s) : s(s) {}
        const SimpleDAQ::TriggerStatus& s;

        typedef Indexed var;
        typedef LogicExpression<Indexed> expr;

        //
        void operator()(var& v) const {
            if (s[v.index]) v.count++;
        }

        void operator()(binop<op_and, Indexed>& b) const { increment(b.oper1, b.oper2); }
        void operator()(binop<op_or, Indexed>& b) const { increment(b.oper1, b.oper2); }
        void operator()(binop<op_xor, Indexed>& b) const { increment(b.oper1, b.oper2); }

        void increment(expr& l, expr& r) const {
            boost::apply_visitor(*this, l);
            boost::apply_visitor(*this, r);
        }

        void operator()(unop<op_not, Indexed>& u) const {
            boost::apply_visitor(*this, u.oper1);
        }
    };

    struct Eval : boost::static_visitor<bool>
    {
        Eval(const SimpleDAQ::TriggerStatus& s) : s(s) {}
        const SimpleDAQ::TriggerStatus& s;

        typedef Indexed var;
        typedef LogicExpression<Indexed> expr;

        //
        bool operator()(const var& v) const {
            return s[v.index] && (v.count % v.reduction == 0);
        }

        bool operator()(const binop<op_and, Indexed>& b) const {
            return boost::apply_visitor(*this, b.oper1) && boost::apply_visitor(*this, b.oper2);
        }
        bool operator()(const binop<op_or, Indexed>& b) const {
            return boost::apply_visitor(*this, b.oper1) || boost::apply_visitor(*this, b.oper2);
        }
        bool operator()(const binop<op_xor, Indexed>& b) const {
            return boost::apply_visitor(*this, b.oper1) ^ boost::apply_visitor(*this, b.oper2);
        }

        bool operator()(const unop<op_not, Indexed>& u) const {
            return !boost::apply_visitor(*this, u.oper1);
        }
    };


    struct Printer : boost::static_visitor<std::string>
    {
        Printer() = default;

        typedef Indexed var;
        typedef LogicExpression<Indexed> expr;

        //
        std::string operator()(const var& v) const {
            return v.name + "(" + to_string(v.reduction) + ")";
        }

        std::string operator()(const binop<op_and, Indexed>& b) const {
            return print(b.oper1, b.oper2, "&");
        }

        std::string operator()(const binop<op_or, Indexed>& b) const {
            return print(b.oper1, b.oper2, "|");
        }

        std::string operator()(const binop<op_xor, Indexed>& b) const {
            return print(b.oper1, b.oper2, "^");
        }

        std::string print(const expr& l, const expr& r, std::string op) const {
            auto left = boost::apply_visitor(*this, l);
            auto right = boost::apply_visitor(*this, r);
            return fmt::format("( {} {} {} )", left, op, right);
        }

        std::string operator()(const unop<op_not, Indexed>& u) const {
            auto inner = boost::apply_visitor(*this, u.oper1);
            return fmt::format("!({})", inner);
        }
    };
}

BOOST_FUSION_ADAPT_STRUCT(
        Indexed,
(std::string, name)
(std::size_t, reduction)
)

namespace simX {
    namespace parser {
        SimpleDAQ::TriggerFunction parseTriggerFunction(const string &s,
                                                          const DetectionSystem &detectionSystem) {
            vector<string> names;
            for (auto& p : detectionSystem.getDetectors()) {
                names.push_back(p->getName());
            }

            return parseTriggerFunction(s, names);
        }

        daq::SimpleDAQ::TriggerFunction parseTriggerFunction(const std::string &s,
                                                             const std::vector<std::string> &names) {
            auto f(std::begin(s)), l(std::end(s));
            DownscaledGrammar<Indexed, decltype(f)> p;

            LogicExpression<Indexed> result;
            bool ok = qi::phrase_parse(f,l,p >> qi::eoi,qi::space,result);
            if (!ok) throw std::invalid_argument("failed to parse '" + s + "'");


            boost::apply_visitor(InitStructures(names), result);

            auto expr = boost::apply_visitor(Printer{}, result);
            auto logger = log::getLogger("DAQTriggerParser");
            logger->debug(expr);

            return [=](const SimpleDAQ::TriggerStatus& status) mutable {
                boost::apply_visitor(IncrementCounts(status), result);
                return boost::apply_visitor(Eval(status), result);
            };
        }
    }
}


