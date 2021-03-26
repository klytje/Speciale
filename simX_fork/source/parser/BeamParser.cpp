//
// Created by munk on 13-12-15.
//

#include "simX/parser/BeamParser.h"
#include "simX/parser/grammar/BeamGrammar.h"

#include <ausa/util/FileUtil.h>
#include <ausa/util/memory>
#include <memory>

using namespace std;
using namespace simX;
using namespace simX::parser;

namespace {
    using DistributionFactories = std::map<std::string, BeamParser::DistributionFactory>;

    unique_ptr<TF1> buildGaussian(BeamParser::Options& opt, AUSA::UnitParser& p) {
        if (!opt.count("sigma"))
            throw std::invalid_argument("Gaussian must have sigma!");

        auto sigma = p.parse(opt["sigma"]);

        if (sigma == 0) return nullptr;

        auto tf1 = make_unique<TF1>("", "gaus(0)", -3 * sigma, 3 * sigma);
        tf1->SetParameters(1, 0, sigma);
        return move(tf1);
    }

    unique_ptr<TF1> buildUniform(BeamParser::Options& opt, AUSA::UnitParser& p) {
        if (!opt.count("low"))
            throw std::invalid_argument("Uniform must have 'low'!");

        if (!opt.count("high"))
            throw std::invalid_argument("Uniform must have 'high'!");

        auto low = p.parse(opt["low"]);
        auto high = p.parse(opt["high"]);

        return (low != high) ? make_unique<TF1>("", "1", low, high) : nullptr;
    }

    unique_ptr<TF1> buildTF1(BeamParser::Options& opt, AUSA::UnitParser& p) {
        if (!opt.count("low"))
            throw std::invalid_argument("TF1 must have 'low'!");

        if (!opt.count("high"))
            throw std::invalid_argument("TF1 must have 'high'!");

        if (!opt.count("expression"))
            throw std::invalid_argument("TF1 must have 'expression'!");

        auto low = p.parse(opt["low"]);
        auto high = p.parse(opt["high"]);

        return (low != high) ? make_unique<TF1>("", opt["expression"].c_str(), low, high) : nullptr;
    }


    unique_ptr<TF1> buildPoint(BeamParser::Options& opt, AUSA::UnitParser& p) {
        return nullptr;
    }


    // XY dist

    unique_ptr<TF1> buildXyGaussian(BeamParser::Options& opt, AUSA::UnitParser& p) {
        auto f = buildGaussian(opt, p);
        if (f) f->SetRange(0, f->GetXmax()); // guard against nullptr
        return move(f);
    }

    unique_ptr<TF1> buildXyUniform(BeamParser::Options& opt, AUSA::UnitParser& p) {
        if (!opt.count("radius"))
            throw std::invalid_argument("XY uniform must have 'radius'!");

        auto radius = std::abs(p.parse(opt["radius"]));
        return (radius > 0) ? make_unique<TF1>("", "1",0.,radius) : nullptr;
    }

    unique_ptr<TF1> buildXyTF1(BeamParser::Options& opt, AUSA::UnitParser& p) {
        if (!opt.count("radius"))
            throw std::invalid_argument("TF1 must have 'radius'!");

        if (!opt.count("expression"))
            throw std::invalid_argument("TF1 must have 'expression'!");

        auto radius = std::abs(p.parse(opt["radius"]));

        return (radius > 0) ? make_unique<TF1>("", opt["expression"].c_str(), 0, radius) : nullptr;
    }

    unique_ptr<TF1> buildDistribution(std::pair<std::string, boost::optional<KeyValueMap>>& val, DistributionFactories& factories, AUSA::UnitParser& p) {
        BeamParser::Options empty;
        auto& name = val.first;
        auto& opt = val.second.get_value_or(empty);

        if (!factories.count(name))
            throw std::invalid_argument("Missing factory for distribution '" + name + "'");


        auto& f = factories[name];
        return f(opt, p);
    }


    Beam buildBeam(BeamPrototype& proto, DistributionFactories& factories, DistributionFactories& xyFactories) {
        if (!proto.energy.is_initialized())
            throw std::invalid_argument("Beam must have an energy");

        if (!proto.center.is_initialized())
            throw std::invalid_argument("Beam must have a center!");

        auto E = proto.energy.get();
        auto theta = proto.theta.get_value_or(0);
        auto phi = proto.phi.get_value_or(0);
        auto center = (proto.center.is_initialized()) ? proto.center.get() : boost::make_tuple(0., 0., 0.);


        unique_ptr<TF1> thetaDist, phiDist, energyDist, xyDist;

        if (proto.thetaDist.is_initialized()) {
            AUSA::UnitParser parser(' ', "deg");
            thetaDist = buildDistribution(proto.thetaDist.get(), factories, parser);
        }

        if (proto.phiDist.is_initialized()) {
            AUSA::UnitParser parser(' ', "deg");
            phiDist = buildDistribution(proto.phiDist.get(), factories, parser);
        }

        if (proto.energyDist.is_initialized()) {
            AUSA::UnitParser parser('k', "eV");
            energyDist = buildDistribution(proto.energyDist.get(), factories, parser);
        }

        if (proto.xyDist.is_initialized()) {
            AUSA::UnitParser parser('m', "m");
            xyDist = buildDistribution(proto.xyDist.get(), xyFactories, parser);
        }


        return Beam(E.first, get<0>(center), get<1>(center), get<2>(center), theta, phi, move(energyDist), move(xyDist), move(thetaDist), move(phiDist), E.second);
    }
}

BeamParser::BeamParser() {
    registerFactory("GAUSSIAN", buildGaussian);
    registerFactory("UNIFORM", buildUniform);
    registerFactory("POINT", buildPoint);
    registerFactory("TF1", buildTF1);

    registerXYFactory("GAUSSIAN", buildXyGaussian);
    registerXYFactory("UNIFORM", buildXyUniform);
    registerXYFactory("POINT", buildPoint);
    registerXYFactory("TF1", buildXyTF1);
}

void BeamParser::registerFactory(std::string key, BeamParser::DistributionFactory f) {
    distFactories[key] = f;
}

void BeamParser::registerXYFactory(std::string key, BeamParser::DistributionFactory f) {
    xyFactories[key] = f;
}

Beam BeamParser::parseFile(const std::string &file) {
    return parse(AUSA::asciiFileToString(file));
}

Beam BeamParser::parse(const std::string &s) {
    using boost::spirit::qi::phrase_parse;
    using boost::spirit::ascii::space;

    BeamPrototype result;
    auto begin = s.begin();
    auto last = s.end();
    bool r = phrase_parse(begin, last, BeamGrammar{}, space, result);
    if (!r) throw std::invalid_argument("failed to parse '" + s + "'");
    if (begin != last) throw std::runtime_error("Boost did not parse the entire string!");

    return buildBeam(result, distFactories, xyFactories);
}
