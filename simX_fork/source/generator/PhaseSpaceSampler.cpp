//
// Created by munk on 6/16/16.
//

#include "simX/generator/PhaseSpaceSampler.h"
#include "simX/samplers/SampleSet.h"
#include "simX/samplers/TextSampleSet.h"
#include "simX/samplers/RootSampleSet.h"
#include "simX/samplers/SampleSampler.h"
#include "simX/samplers/MultiSampler.h"

#include <ausa/util/UnitParser.h>

#include <cmath>
#include <ausa/util/FileUtil.h>
#include <ausa/util/StringUtil.h>
#include <simX/samplers/JonasLibSampleSet.h>
#include <simX/Logger.h>

using namespace simX;

PhaseSpaceSampler::PhaseSpaceSampler(Input input_, std::vector<Particle *> daughters, std::string unit)
    : daughters(daughters), result(daughters.size()), input(move(input_))
{
    if (!input->hasNext()) std::invalid_argument("SampleSet to PhaseSpaceSampler must not be empty!");
    if (input->isSampled()) std::invalid_argument("SampleSet to PhaseSpaceSampler is not sampled!");

    AUSA::UnitParser parser('k', "eV");
    multiplier = parser.parse(std::to_string(1) + unit);
}

std::vector<TLorentzVector> &PhaseSpaceSampler::getFourMomenta() {
    static auto LOGGER = simX::log::getLogger("PhaseSpaceSampler");

    if (!input->hasNext()) input->reset();
    input->next();
    auto& sample = input->getSample();

    for (size_t i = 0; i < daughters.size(); i++) {
        auto& p = result[i];
        auto  d = daughters[i];
        p.SetE(0);
        p.SetPx(multiplier*sample[i*3+0]);
        p.SetPy(multiplier*sample[i*3+1]);
        p.SetPz(multiplier*sample[i*3+2]);

        double mx = d->getMass() + d->getExcitationEnergy();
        p.SetE(std::sqrt(std::pow(mx, 2) - p.M2()));
        LOGGER->trace("{} (E,px,py,pz) = ({},{},{},{})", d->getName(), p.E(), p.X(), p.Y(), p.Z());
    }

    return result;
}

std::unique_ptr<FinalStateGenerator>
PhaseSpaceSampler::build(const NBodyDecay &d, parser::ReactionParser::Options &opt) {
    if (!opt.count("file")) throw std::invalid_argument("File must be specified for PhaseSpaceSampler");
    if (!opt.count("unit")) throw std::invalid_argument("Unit must be specified for PhaseSpaceSampler");

    auto& file = opt["file"];
    auto& unit = opt["unit"];
    auto n = d.getDaughters().size();


    auto files = AUSA::findFilesMatchingWildcard(file);
    if (files.empty())
        throw std::invalid_argument("Could not find files matching " + file);

    std::vector<std::unique_ptr<samplers::SampleSet>> generators;
    for (auto& single : files) {
        std::unique_ptr<samplers::SampleSet> generator;
        if (AUSA::endswith(single, ".root")) {

            TDirectory::TContext{nullptr};
            TFile f{single.c_str()};
            if (f.Get("sim"))
                generator = std::make_unique<samplers::JonasLibSampleSet>(single, n);
            else if (f.Get("sample"))
                generator = std::make_unique<samplers::RootSampleSet>(single, 3*n);
            else
                throw std::invalid_argument("Unknown ROOT based PS file");

        } else {
            generator = std::make_unique<samplers::TextSampleSet>(single, 3*n);
        }

        generators.push_back(move(generator));
    }

    std::unique_ptr<samplers::SampleSet> generator;
    if (generators.size() == 1) generator = move(generators.front());
    else generator = std::make_unique<samplers::MultiSampler>(move(generators));

    if (!generator->isSampled()) generator = std::make_unique<samplers::SampleSampler>(move(generator));

    return std::make_unique<PhaseSpaceSampler>(move(generator), d.getDaughters(), unit);
}
