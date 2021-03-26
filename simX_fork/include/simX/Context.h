//
// Created by munk on 07-11-17.
//

#ifndef SIMX_CONTEXT_H
#define SIMX_CONTEXT_H

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
#define SIMX_TCONTEXT TDirectory::TContext __________simx_context;
#else
#define SIMX_TCONTEXT TDirectory::TContext __________simx_context{nullptr};
#endif


#endif //SIMX_CONTEXT_H
