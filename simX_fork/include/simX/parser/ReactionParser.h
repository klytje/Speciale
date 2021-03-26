//
// Created by munk on 28-06-15.
//

#ifndef SIMX_REACTIONPARSER_H
#define SIMX_REACTIONPARSER_H

#include "simX/ProcessChain.h"
#include "simX/NBodyDecay.h"

#include <string>
#include <memory>
#include <map>
#include <functional>


namespace simX {
    namespace angular {
        class AngularCorrelation;
    }


    namespace parser {

        /**
         * Object responsible for parsing a reaction file to a ProcessChain.
         */
        class ReactionParser {
        public:

            using ChainPtr = std::unique_ptr<ProcessChain>;

            struct ReactionBranch {
                ChainPtr chain;
                double ratio;
                std::string name;
                std::string path;
            };




            /**
             * Lower and upper limits.
             */
            using Limits = std::pair<double, double>;

            /**
             * A string -> string map with additional user options.
             */
            using Options = std::map<std::string, std::string>;

            /**
             * A function object which will construct an instance of AngularCorrelation.
             * @param phi The lower and upper limit of the phi area the user want to sample.
             * @param theta The lower and upper limit of the theta area the user want to sample.
             * @param opt A map with additional options the user have specified.
             */
            using AngularFactory = std::function<std::unique_ptr<angular::AngularCorrelation>(const Limits& phi, const Limits& theta, Options& opt)>;

            using GeneratorFactory = std::function<std::unique_ptr<FinalStateGenerator>(const NBodyDecay& d, Options& opt)>;

            using WeightFactory = std::function<std::unique_ptr<WeightCalculator>(NBodyDecay& d, Options& opt)>;

            /**
             * Construct the default ReactionParser.
             */
            ReactionParser();
            ~ReactionParser();

            /**
             * Parse the supplied input and return a unique_ptr to a ProcessChain.
             *
             * @return A unique_ptr to a ProcessChain.
             * @throw std::invalid_argument if there is mistakes in the supplied input.
             */
            ChainPtr parseString(const std::string& input);

            /**
             * Open the file specified by the string and parse the file.
             *
             * @copydoc parseString
             */
            ChainPtr parseFile(const std::string& input);


            std::vector<ReactionBranch> parseMultiple(const std::string& file);

            /**
             * Register a AngularFactory with a identifier.
             * If the user supplies this exact identifier then your AngularFactory
             * is asked to construct the AngularCorrelation object.
             *
             * @param identifier The identifier for which you want to be called.
             * @param f The factory object which will construct your AngularCorrelation.
             */
            void registerFactory(const std::string& identifier, AngularFactory f);

            /**
             * Register a GeneratorFactory with a identifier.
             * If the user supplies this exact identifier then your GeneratorFactory
             * is asked to construct the FinalStateGenerator object.
             *
             * @param identifier The identifier for which you want to be called.
             * @param f The factory object which will construct your FinalStateGenerator.
             */
            void registerFactory(const std::string& identifier, GeneratorFactory f);

            void registerFactory(const std::string& identifier, WeightFactory f);

        private:
            // We use Pimpl not to expose ReactionGrammar which takes forever to compile.
            struct Pimpl;
            std::unique_ptr<Pimpl> pimpl;
        };
    }
}
#endif //SIMX_REACTIONPARSER_H
