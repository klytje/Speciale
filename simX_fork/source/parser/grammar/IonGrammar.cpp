//
// Created by munk on 25-06-15.
//

#include "simX/parser/grammar/IonGrammar.h"

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix.hpp>

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

simX::parser::ElementSymbolTable::ElementSymbolTable() {
    add
            ("H"	,	1)
            ("He"	,	2)
            ("Li"	,	3)
            ("Be"	,	4)
            ("B"	,	5)
            ("C"	,	6)
            ("N"	,	7)
            ("O"	,	8)
            ("F"	,	9)
            ("Ne"	,	10)
            ("Na"	,	11)
            ("Mg"	,	12)
            ("Al"	,	13)
            ("Si"	,	14)
            ("P"	,	15)
            ("S"	,	16)
            ("Cl"	,	17)
            ("Ar"	,	18)
            ("K"	,	19)
            ("Ca"	,	20)
            ("Sc"	,	21)
            ("Ti"	,	22)
            ("V"	,	23)
            ("Cr"	,	24)
            ("Mn"	,	25)
            ("Fe"	,	26)
            ("Co"	,	27)
            ("Ni"	,	28)
            ("Cu"	,	29)
            ("Zn"	,	30)
            ("Ga"	,	31)
            ("Ge"	,	32)
            ("As"	,	33)
            ("Se"	,	34)
            ("Br"	,	35)
            ("Kr"	,	36)
            ("Rb"	,	37)
            ("Sr"	,	38)
            ("Y"	,	39)
            ("Zr"	,	40)
            ("Nb"	,	41)
            ("Mo"	,	42)
            ("Tc"	,	43)
            ("Ru"	,	44)
            ("Rh"	,	45)
            ("Pd"	,	46)
            ("Ag"	,	47)
            ("Cd"	,	48)
            ("In"	,	49)
            ("Sn"	,	50)
            ("Sb"	,	51)
            ("Te"	,	52)
            ("I"	,	53)
            ("Xe"	,	54)
            ("Cs"	,	55)
            ("Ba"	,	56)
            ("La"	,	57)
            ("Ce"	,	58)
            ("Pr"	,	59)
            ("Nd"	,	60)
            ("Pm"	,	61)
            ("Sm"	,	62)
            ("Eu"	,	63)
            ("Gd"	,	64)
            ("Tb"	,	65)
            ("Dy"	,	66)
            ("Ho"	,	67)
            ("Er"	,	68)
            ("Tm"	,	69)
            ("Yb"	,	70)
            ("Lu"	,	71)
            ("Hf"	,	72)
            ("Ta"	,	73)
            ("W"	,	74)
            ("Re"	,	75)
            ("Os"	,	76)
            ("Ir"	,	77)
            ("Pt"	,	78)
            ("Au"	,	79)
            ("Hg"	,	80)
            ("Tl"	,	81)
            ("Pb"	,	82)
            ("Bi"	,	83)
            ("Po"	,	84)
            ("At"	,	85)
            ("Rn"	,	86)
            ("Fr"	,	87)
            ("Ra"	,	88)
            ("Ac"	,	89)
            ("Th"	,	90)
            ("Pa"	,	91)
            ("U"	,	92)
            ("Np"	,	93)
            ("Pu"	,	94)
            ("Am"	,	95)
            ("Cm"	,	96)
            ("Bk"	,	97)
            ("Cf"	,	98)
            ("Es"	,	99)
            ("Fm"	,	100)
            ("Md"	,	101)
            ("No"	,	102)
            ("Lr"	,	103)
            ("Rf"	,	104)
            ("Db"	,	105)
            ("Sg"	,	106)
            ("Bh"	,	107)
            ("Hs"	,	108)
            ("Mt"	,	109)
            ("Ds"	,	110)
            ("Rg"	,	111)
            ("Cp"	,	112)
            ("Uut"	,	113)
            ("Uuq"	,	114)
            ("Uup"	,	115)
            ("Uuh"	,	116)
            ("Uus"	,	117)
            ("Uuo"	,	118)
            ;
}

simX::parser::AzeIsotopeGrammar::AzeIsotopeGrammar()  : AzeIsotopeGrammar::base_type(start) {
    using qi::eps;
    using qi::lit;
    using qi::_val;
    using qi::uint_;
    using qi::_1;
    using AUSA::EnergyLoss::Ion;



    start = (element_name >> -lit('-') >> uint_)
            [ qi::_val = boost::phoenix::construct<Ion>(qi::_1, qi::_2) ]
            |
            (uint_ >> element_name)
            [ qi::_val = boost::phoenix::construct<Ion>(qi::_2, qi::_1) ]
            ;

}

simX::parser::PredefinedIons::PredefinedIons() {
    using AUSA::EnergyLoss::Ion;
    add
            ("g"	,	Ion(0,0))
            ("n"	,	Ion(0,1))
            ("p"	,	Ion(1,1))
            ("d"	,	Ion(1,2))
            ("t"	,	Ion(1,3))
            ("a"	,	Ion(2,4))
            ;
}

simX::parser::IonGrammar::IonGrammar() : IonGrammar::base_type(start) {
    start %= aze | predefined;
}
