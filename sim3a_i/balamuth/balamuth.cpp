#include "cml.h"
#include "Decay.h"
#include "TextWriter.h"
#include "RootInterferenceWriter.h"

#include <iostream>
#include <memory>
#include <ausa/util/memory>
#include <ausa/util/FileUtil.h>
#include <ausa/util/StringUtil.h>
#include <ausa/parser/UnitParser.h>
#include <string>

#include "ausa/constants/Mass.h"
#include "jonas/TripleDecay.h"
#include "jonas/Pair.h"
#include "jonas/Nucleus.h"
#include "jonas/Level.h"
#include "jonas/Channel.h"
#include "jonas/BalamuthWeight.h"
#include "jonas/BalamuthInterferenceWeight.h"
#include "jonas/TripleInterferenceDecay.h"

using namespace std;
using namespace AUSA::Constants;

using FillFunction = std::function<void(const TriAlphaDecay&)>;

unique_ptr<BalamuthWeight> createBalamuthWeight(CommandLineOptions& cml,int L) {
    auto balamuth = std::make_unique<BalamuthWeight>();

    Nucleus He4(4,2,0,1);
    Nucleus Be8(8,4,0,1);

    // First pair
    Pair pPair(Be8,He4);

    // Carbon system
    auto pChannel = std::make_unique<Channel>(pPair,L,35. /*Arb*/ ,cml.r0p, true);
    Level pLevel(42 /*Arb*/ ,cml.J1,1);

    //The second channel
    Pair sPair(He4,He4);
    auto sChannel = std::make_unique<Channel>(sPair,cml.J2,cml.gs,cml.r0s, true);
    Level sLevel(cml.ExBe,cml.J2,1);

	std::cout << sLevel.E() << "\n";
	double gamma = pow(36,2);
	std::cout << gamma/(1+gamma*sChannel->ShiftFunctionDeriv(sLevel.E())) << std::endl;
	std::cout << 2*gamma/(1+gamma*sChannel->ShiftFunctionDeriv(sLevel.E()))*sChannel->Penetrability(sLevel.E()) << std::endl;
	std::flush(std::cout);
	double Er = sLevel.E();
	for (double Er = sLevel.E()-500; Er<sLevel.E()+501; Er+=10) {
	//	std::cout << Er << " " << sLevel.E()-pow(sChannel->ReducedWidth(),2)*sChannel->Penetrability(Er)-Er << " " << sChannel->Penetrability(Er) << "\n";
	}

    balamuth->SetPrimary(pLevel, move(pChannel));
    balamuth->SetSecondary(sLevel, move(sChannel));

    balamuth->DoCorrection(cml.coloumbCorrection, cml.coloumbRadius);
    return move(balamuth);
}

void generateDecays(const double Q, const FillFunction& fill, CommandLineOptions& cml) {
    TripleInterferenceDecay decay;
    decay.SetQ(Q);
    decay.DoRecoil(false);
    if (cml.ExBe < 5) decay.DoGroundState(true);

    if (!cml.interference) {
        auto balamuth = createBalamuthWeight(cml,cml.L);
        decay.SetDecayWeight(move(createBalamuthWeight(cml,cml.L)));
    } else {
        auto balamuth1 = createBalamuthWeight(cml,cml.L);
        auto balamuth2 = createBalamuthWeight(cml,cml.L2);
        auto balamuth = new BalamuthInterferenceWeight(move(balamuth1),move(balamuth2));
        decay.SetDecayWeights({balamuth});
    }

    auto printFrac = max(cml.N/1000,1UL);

    TriAlphaDecay info{};

    size_t discard = 0;

    for (size_t i = 0; i < cml.N; ++i) {
        if (i % printFrac == 0) cout << i << "\r" << flush;

        info.w = decay.Generate();

        info.wB = decay.GetCalculatorWeight();
        info.wU = decay.GetGeneratorWeight();
        info.factors = decay.GetFactors();

        for (int j = 0; j < 3; ++j) {
            info.eCM[j] = decay.GetCmEnergy(j);
        }
        TMath::Sort(3,info.eCM,info.order,true);

        for (int j = 0; j < 3; ++j) {
            info.p[j] = &decay.GetProduct(info.order[j]);
        }


        if (std::isnan(info.factors[0][0])) cout << "ERR!!!\n";

        if (std::isinf(info.w) || std::isnan(info.w) || std::isnan(info.factors[0][0]+info.factors[0][1]+info.factors[0][2]+info.factors[0][3]) || std::isinf(info.factors[0][0]+info.factors[0][1]+info.factors[0][2]+info.factors[0][3])) {
            discard++;
            cerr << "Invalid weight. Discarding event!" << endl;
            i--;

            if (discard > cml.N/1000.) {
                cerr << "Weight of more than 1:1000 is invalid. Stopping!" << endl;
                exit(255);
            }
        }
        else {
            fill(info);
        }
    }
    cout << cml.N << "\r" << endl;
}


int main(int argc, char* argv[]) {
    auto cml = parseCml(argc, argv);

    cout << "C12 excitation " << cml.ExC << " keV" << endl;
    auto Q = (isotopeMass("C12") + cml.ExC - 3 * isotopeMass("He4"));

    cout << "Q = " << Q << " keV" << endl;

    std::unique_ptr<TextWriter> text;
    std::unique_ptr<RootInterferenceWriter> root;
    FillFunction fill;

    const auto textOutput = !AUSA::endswith(cml.output, ".root");
    if (textOutput) {
        text = std::make_unique<TextWriter>(cml.output);
        fill = [&](const TriAlphaDecay& decay) {text->fill(decay);};
    }
    else {
        root = std::make_unique<RootInterferenceWriter>(cml.output);
        fill = [&](const TriAlphaDecay& decay) {root->fill(decay);};
    }

    generateDecays(Q, fill,cml);

}