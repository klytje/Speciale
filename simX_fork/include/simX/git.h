//
// Created by munk on 31-08-15.
//

#ifndef SIMX_GIT_H
#define SIMX_GIT_H
namespace simX {
    /**
     * Name of the checked out branch used when compiling simX.
     */
    extern const char GIT_BRANCH[];

    /**
     * The sha-1 value associated with the latest commit on simX.
     */
    extern const char GIT_HASH[];
}
#endif //SIMX_GIT_H
