//
// Created by munk on 07-05-15.
//


#include <unittest++/UnitTest++.h>
#include <unittest++/TestRunner.h>
#include <unittest++/TestReporterStdout.h>
#include <string.h>
#include <unordered_set>

using namespace UnitTest;

int main(int argc, const char * argv[])
{
    if (argc > 1 && strcmp(argv[1], "--all") == 0)
        return UnitTest::RunAllTests();

    std::unordered_set<std::string> disabled = {
            "MCStragglingTest"
    };

    TestReporterStdout reporter;
    TestRunner runner(reporter);
    return runner.RunTestsIf(Test::GetTestList(), nullptr,
                             [&](Test* t) {
                                 return !disabled.count(t->m_details.suiteName);
                                 return strcmp(t->m_details.suiteName, "MultiSamplerTest") == 0;
//                                 return strcmp(t->m_details.testName, "Outside36DegreesStrip3GivesNothing") == 0;
                             }, 0);
}
