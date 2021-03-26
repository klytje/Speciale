//
// Created by munk on 28-06-15.
//

#include "simX/parser/ReactionParser.h"
#include "simX/parser/grammar/ReactionGrammar.h"
#include "simX/CompoundFormation.h"
#include "simX/TwoBodyDecay.h"

#include "simX/angular/AngularCorrelationTF1.h"
#include "simX/angular/IsotropicAngularCorrelation.h"
#include "simX/angular/FixedAngularCorrelation.h"
#include "simX/angular/GraphDistribution.h"
#include "simX/angular/RutherfordDistribution.h"

#include <ausa/util/memory>
#include <ausa/util/FileUtil.h>

#include <iostream>
#include <unordered_map>
#include <simX/generator/NBodyPhaseSpace.h>
#include <simX/generator/TwoBodyKinematics.h>
#include <simX/generator/SingleBodyKinematics.h>
#include "simX/generator/PhaseSpaceSampler.h"
#include <simX/weight/GammaPhaseSpace.h>
#include <simX/weight/BreitWignerWithPenetrability.h>
#include <simX/weight/UniformWeightCalculator.h>
#include <simX/weight/SimpleBreitWigner.h>
#include <simX/weight/ProductWeightCalculator.h>
#include <simX/weight/TF1WeightCalculator.h>

#include <boost/algorithm/string.hpp>
#include <ausa/json/MalformedJsonException.h>
#include <rapidjson/document.h>
#include <ausa/json/JSONUtil.h>
#include "ausa/util/FileUtil.h"

using namespace std;
using namespace simX;
using namespace simX::parser;
using namespace simX::angular;

namespace {
    simX::Particle convert(const simX::parser::Particle& p, simX::Particle* parent) {
        return simX::Particle{p.ion, p.doTracking(), p.ex(), p.G0(), parent};
    }
}

struct ReactionParser::Pimpl {
    using OptReaction = boost::optional<ReactionDescription>;
    using PhysParticle = simX::Particle;
    using DecayPtr = ProcessChain::DecayPtr;
    using AdPtr = std::unique_ptr<AngularCorrelation>;
    using GeneratorPtr = unique_ptr<FinalStateGenerator>;
    using WPtr = unique_ptr<WeightCalculator>;



    std::vector<ReactionBranch> parseMultiple(const std::string& file) {
        using namespace rapidjson;
        using namespace AUSA::JSON;
        using namespace AUSA;
        std::vector<ReactionBranch> chains;
        std::string input = AUSA::asciiFileToString(file);

        Document document;
        // If we can't parse with json, we assume it is a normal reaction file and
        // process as such
        try {
            tryParse(input, document);
        } catch (MalformedJsonException e) {
            chains.push_back({ReactionParser().parseFile(file), 1.0, file, file});
            return chains;
        }

        // We read in all the files and normalize their ratios
        double totalRation = 0.0;
        for(auto i=document.MemberBegin(); i< document.MemberEnd(); ++i) {
            double ratioValue = i->value.GetDouble();
            totalRation += ratioValue;
            chains.push_back({
                                     ReactionParser().parseFile(i->name.GetString()),
                                     ratioValue,
                                     extractFileName(i->name.GetString()),
                                     i->name.GetString()
                             });
        }
        // Normalize ratios
        for(ReactionBranch& b : chains) {
            b.ratio /= totalRation;
        }
        return chains;
    }

    ChainPtr parse(const std::string& input) {
        using boost::spirit::qi::phrase_parse;
        using boost::spirit::ascii::space;


        ReactionChainDescription result;
        bool r = phrase_parse(begin(input), end(input), ReactionChainGrammar{}, space, result);
        if (!r) throw std::invalid_argument("Failed to parse input !");

        ProcessChain::Tree t;
        auto top = t.begin();
        
        auto beam = make_unique<PhysParticle>(convert(result.beam, nullptr));
        auto first = beam.get();

        if (result.target) {
            auto fusion = make_unique<CompoundFormation>(*beam, convert(result.target.get(), nullptr));
            first = fusion -> getDaughters()[0];
            top = t.insert(top, std::move(fusion));
        }

        parseReaction(result.chain, first, top, t);

        return make_unique<ProcessChain>(move(beam), t);
    }

    void parseReaction(OptReaction &r, PhysParticle* parent, ProcessChain::Tree::iterator iter, ProcessChain::Tree &tree) {
        if (!r) return;

        ReactionDescription& reaction = r.get();

        auto& prod = reaction.products;

        if (prod.size() < 1)
            throw invalid_argument("NuclearProcess of " + parent -> getName() + " to " + to_string(prod.size()) + " daughters is not allowed!");

        // Make decay object
        DecayPtr decay = makeDecay(reaction, parent);

        // Insert into tree
        if (tree.empty())
            iter = tree.insert(iter, std::move(decay));
        else
            iter = tree.append_child(iter, std::move(decay));

        // Now we make all the daughters decay
        for (size_t i = 0; i < prod.size(); i++) {
            auto p = (*iter) -> getDaughters()[i];
            parseReaction(prod[i].second, p, iter, tree);
        }
    }

    DecayPtr makeDecay(ReactionDescription& reaction, PhysParticle* parent) {
        auto& s = reaction.settings;

        vector<simX::Particle> daughters;
        for (auto& d : reaction.products) daughters.emplace_back(convert(d.first, parent));

        auto ac = makeAngularCorrelation(reaction);
        auto L = (s.is_initialized()) ? s -> L.get_value_or(0) : 0;

        auto decay = make_unique<NBodyDecay>(*parent, daughters, nullptr, move(ac), nullptr, L);

        auto generator = makeGenerator(*decay, *parent, reaction);
        auto weight = makeWeight(*decay, *parent, reaction);

        decay -> setFinalStateGenerator(move(generator));
        decay -> setWeightCalculator(move(weight));

        return move(decay);
    }

    WPtr makeWeight(NBodyDecay& np, PhysParticle& parent, ReactionDescription& r) {
        auto& daughters = np.getDaughters();
        auto n = daughters.size();

        bool gamma = false;
        for (auto d : daughters) gamma |= (d->getName() == "gamma");

        auto& s = r.settings;
        if (!(s.is_initialized() && s->weight.is_initialized())) { // We apply defaults!
            if (n == 2) {

                auto bw = make_unique<BreitWignerWithPenetrability>(np);

                if (gamma)
                    return make_unique<ProductWeightCalculator>(move(bw), make_unique<GammaPhaseSpace>(np));
                else
                    return move(bw);
            }

            return make_unique<SimpleBreitWigner>(np);
        }

        auto& wId = s->weight->first;

        std::vector<std::string> ids;
        boost::split(ids, wId, boost::is_any_of("*"));

        WPtr result = makeCustomWeight(ids[0], s->weight->second, np);
        for (size_t i = 1; i < ids.size(); i++) {
            result = make_unique<ProductWeightCalculator>(move(result), makeCustomWeight(ids[i], s->weight->second, np));
        }
        return result;
    }

    WPtr makeCustomWeight(const std::string& id, boost::optional<KeyValueMap>& opt, NBodyDecay& np) {
        if (!weightFactory.count(id))
            throw std::invalid_argument("Unknown weight function '" + id + "' !");

        return weightFactory[id](np, opt.get_value_or(emptyOptions));
    }

    GeneratorPtr makeGenerator(NBodyDecay& np, PhysParticle& parent, ReactionDescription& r) {
        auto& s = r.settings;
        if (s.is_initialized()) {
            if (s->gen.is_initialized()) {
                auto& name = s->gen->first;
                if (!generatorFactory.count(name)) throw runtime_error("Unknown phase space generator: " + name);

                auto& f = generatorFactory[name];
                return f(np, s->gen->second.get_value_or(emptyOptions));
            }
        }

        auto& d = np.getDaughters();
        if (d.size() == 1)
            return make_unique<SingleBodyKinematics>(parent);
        else if (d.size() == 2)
            return make_unique<TwoBodyKinematics>(parent, *d[0], *d[1]);
        else
            return make_unique<NBodyPhaseSpace>(parent, d);
    }

    AdPtr makeAngularCorrelation(ReactionDescription& r) {
        auto& s = r.settings;
        string ad = "ISO";
        auto phi = make_pair<double, double>(0., 360.);
        auto theta = make_pair<double, double>(0.,180.);
        auto& adOpt = emptyOptions;
        if (s.is_initialized()) {
            if (s -> ad)    ad = s -> ad.get().first;
            if (s -> phi)   phi = s -> phi.get();
            if (s -> theta) theta = s -> theta.get();
            if (s -> ad && s -> ad -> second) adOpt = s -> ad -> second.get();
        }


        if (!adFactory.count(ad))
            throw invalid_argument("I do not know AngularCorrelation '" + ad + "'");

        return adFactory[ad](phi, theta, adOpt);
    }

    void registerFactory(const std::string& identifier, AngularFactory f) {
        adFactory[identifier] = f;
    }

    void registerFactory(const std::string& identifier, GeneratorFactory f) {
        generatorFactory[identifier] = f;
    }

    void registerFactory(const std::string& identifier, WeightFactory f) {
        weightFactory[identifier] = f;
    }

private:
    unordered_map<std::string, AngularFactory> adFactory;
    unordered_map<std::string, GeneratorFactory> generatorFactory;
    unordered_map<std::string, WeightFactory> weightFactory;
    Options emptyOptions;
};


ReactionParser::ReactionParser() : pimpl(std::make_unique<Pimpl>()){
    registerFactory("ISO", [](const Limits& phi, const Limits& theta, Options& m) {
        return make_unique<IsotropicAngularCorrelation>(phi, theta);
    });

    registerFactory("RUTHERFORD", [](const Limits& phi, const Limits& theta, Options& m) {
        return make_unique<RutherfordDistribution>(phi, theta);
    });

    registerFactory("GRAPH", [](const Limits& phi, const Limits& theta, Options& m) {
        if (!m.count("file")) throw invalid_argument("You forgot to specify file for GRAPH!");
        string format = "%lg %lg";
        if (m.count("format")) format = m["format"];

        auto dist = make_unique<GraphDistribution>(phi.first, phi.second, theta.first, theta.second, m["file"], format);
        if (m.count("debug")) dist->printData();

        return move(dist);
    });

    registerFactory("TF1", [](const Limits& phi, const Limits& theta, Options& m) {
        if (!m.count("phi")) throw invalid_argument("You forgot to specify phi for TF1");
        if (!m.count("theta")) throw invalid_argument("You forgot to specify theta for TF1");

        return make_unique<AngularCorrelationTF1>(m["theta"], m["phi"], phi.first, phi.second, theta.first, theta.second);
    });

    registerFactory("FIXED", [&](const ReactionParser::Limits& phi, const ReactionParser::Limits& theta, ReactionParser::Options& m) {
        double th = stof(m["theta"]);
        double ph = stof(m["phi"]);
        return std::make_unique<FixedAngularCorrelation>(th,ph);
    });



    // Weights
    registerFactory("BW", [](NBodyDecay& d, Options& opt) {
        return make_unique<SimpleBreitWigner>(d);
    });

    registerFactory("GAMMA", [](NBodyDecay& d, Options& opt) {
        return make_unique<GammaPhaseSpace>(d);
    });

    registerFactory("PEN", [](NBodyDecay& d, Options& opt) {
        return make_unique<BreitWignerWithPenetrability>(d);
    });

    registerFactory("UNIFORM", [](NBodyDecay& d, Options& opt) {
        return make_unique<UniformWeightCalculator>();
    });

    registerFactory("PREGEN", PhaseSpaceSampler::build);

    // Weights
    registerFactory("TF1", [](NBodyDecay& d, Options& opt) {
        return make_unique<TF1WeightCalculator>(d, opt["f"]);
    });
}

ReactionParser::~ReactionParser() {
    // Do nothing, but here we know sizeof(Pimpl)
}

ReactionParser::ChainPtr ReactionParser::parseFile(const std::string& input) {
    return parseString(AUSA::asciiFileToString(input));
}

ReactionParser::ChainPtr ReactionParser::parseString(const std::string& input) {
    return pimpl -> parse(input);
}

std::vector<ReactionParser::ReactionBranch> ReactionParser::parseMultiple(const std::string& input) {
    return pimpl -> parseMultiple(input);
}


void ReactionParser::registerFactory(const std::string& identifier, ReactionParser::AngularFactory f) {
    pimpl -> registerFactory(identifier, move(f));
}

void ReactionParser::registerFactory(const std::string& identifier, ReactionParser::GeneratorFactory f) {
    pimpl -> registerFactory(identifier, move(f));
}

void ReactionParser::registerFactory(const std::string& identifier, ReactionParser::WeightFactory f) {
    pimpl -> registerFactory(identifier, move(f));
}


