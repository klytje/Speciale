//
// Created by munk on 09-11-16.
//

#include "simX/parser/TdcThresholdParser.h"
#include "simX/Context.h"
#include <ausa/util/FileUtil.h>
#include <ausa/json/JSONUtil.h>
#include <iostream>

using namespace rapidjson;
using namespace std;
using namespace simX;
using namespace simX::daq;
using namespace simX::parser;

TdcThresholds simX::parser::parseTdcThresholds(const std::string &input, size_t required) {
    using namespace AUSA::JSON;

    Document document;
    tryParse(input, document);

    SIMX_TCONTEXT;
    auto rootFile = readString(document, "file");
    TFile file(rootFile.c_str());

    // Stop ROOT for trying to GC these histograms...

    if (!file.IsOpen()) throw std::invalid_argument("Failed opening root file: " + rootFile);

    vector<TH1D> histograms;
    auto& histInput = readArray(document, "histograms");

    if (histInput.Size() != required)
        throw std::invalid_argument("Require " + to_string(required) +
                                            " tdc histograms. Only " + to_string(histInput.Size()) + " provided");

    for (SizeType i = 0; i < histInput.Size(); ++i) {
        auto& entry = histInput[i];
        if (!entry.IsString())
            throw std::invalid_argument("Entry " + to_string(i) + " of array 'histograms is not a string!'");

        auto name = string(entry.GetString());
        if (name == "-DISABLED-") histograms.emplace_back("", "", 1, 0, 100E20);
        else {
            auto h = dynamic_cast<const TH1*>(file.Get(name.c_str()));
            if (!h) throw std::invalid_argument(string("No histogram named '") + name + "' in " + file.GetName());

            histograms.emplace_back("", "", h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            histograms.back().SetDirectory(nullptr);
            for (int j = 1; j <= h->GetNbinsX(); ++j) {
                histograms.back().SetBinContent(i, h->GetBinContent(i));
            }
        }
    }
    file.Close();

    return TdcThresholds(histograms);
}

TdcThresholds simX::parser::parseTdcThresholdsFromFile(const std::string &file, size_t required) {
    try {
        auto str = AUSA::asciiFileToString(file);
        return parseTdcThresholds(str, required);
    } catch (...) {
        cerr << "Failed to parse TDC threshold file " << file << endl;
        throw;
    }
}
