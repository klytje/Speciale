#include "cml.h"
#include "Decay.h"
#include "TextWriter.h"
#include "RootWriter.h"

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

using namespace std;
using namespace AUSA::Constants;

using FillFunction = std::function<void(const TriAlphaDecay&)>;

void generateDecays(const double Q, const size_t N, const int L, const int J1, const int J2, const double r0p, const double r0s,
                    const double gs, const double Es, const bool coloumbCorrection, const double coloumbRadius,
                    const FillFunction& fill) {
    TripleDecay decay;
    decay.SetQ(Q);
    decay.DoRecoil(false);
    if (Es < 5) decay.DoGroundState(true);

    auto balamuth = std::make_unique<BalamuthWeight>();

    Nucleus He4(4,2,0,1);
    Nucleus Be8(8,4,0,1);

    // First pair
    Pair pPair(Be8,He4);

    // Carbon system
    auto pChannel = std::make_unique<Channel>(pPair,L,35. /*Arb*/ ,r0p, true);
    Level pLevel(42 /*Arb*/ ,J1,1);

    //The second channel
    Pair sPair(He4,He4);
    auto sChannel = std::make_unique<Channel>(sPair,J2,gs,r0s, true);
    Level sLevel(Es,J2,1);

    balamuth->SetPrimary(pLevel, move(pChannel));
    balamuth->SetSecondary(sLevel, move(sChannel));

    balamuth->DoCorrection(coloumbCorrection, coloumbRadius);

    decay.SetDecayWeight(move(balamuth));

    auto printFrac = max(N/1000,1UL);

    TriAlphaDecay info{};

    size_t discard = 0;

    for (size_t i = 0; i < N; ++i) {
        if (i % printFrac == 0) cout << i << "\r" << flush;

        info.w = decay.Generate();

        info.wB = decay.GetCalculatorWeight();
        info.wU = decay.GetGeneratorWeight();

        for (int j = 0; j < 3; ++j) {
            info.eCM[j] = decay.GetCmEnergy(j);
        }
        TMath::Sort(3,info.eCM,info.order,true);

        for (int j = 0; j < 3; ++j) {
            info.p[j] = &decay.GetProduct(info.order[j]);
        }

        if (std::isinf(info.w) || std::isnan(info.w)) {
            discard++;
            cerr << "Invalid weight. Discarding event!" << endl;

            if (discard > N/1000.) {
                cerr << "Weight of more than 1:1000 is invalid. Stopping!" << endl;
                exit(255);
            }
        }
        else {
            fill(info);
        }
    }
    cout << N << "\r" << endl;
}


int main(int argc, char* argv[]) {
    auto cml = parseCml(argc, argv);

    cout << "C12 excitation " << cml.ExC << " keV" << endl;
    auto Q = (isotopeMass("C12") + cml.ExC - 3 * isotopeMass("He4"));

    cout << "Q = " << Q << " keV" << endl;

    std::unique_ptr<TextWriter> text;
    std::unique_ptr<RootWriter> root;
    FillFunction fill;

    const auto textOutput = !AUSA::endswith(cml.output, ".root");
    if (textOutput) {
        text = std::make_unique<TextWriter>(cml.output);
        fill = [&](const TriAlphaDecay& decay) {text->fill(decay);};
    }
    else {
        root = std::make_unique<RootWriter>(cml.output);
        fill = [&](const TriAlphaDecay& decay) {root->fill(decay);};
    }



    generateDecays(Q, cml.N, cml.L, cml.J1, cml.J2, cml.r0p, cml.r0s,
                   cml.gs, cml.ExBe, cml.coloumbCorrection, cml.coloumbRadius, fill);
}