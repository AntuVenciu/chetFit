#include <iostream>
#include <unistd.h>
#include <cstdio>

#include "fitteralgorithms.hh"
#include "auxiliaryalgorithms.hh"
#include "options.hh"

using namespace std;


enum class FitterType
{
    PlanarFitter,
    SpacepointFitter,
    HelixFitter,
    Unknown
};



FitterType GetFitterTypeFromString(const std::string& str)
{
    if(str == "planar") return FitterType::PlanarFitter;
    if(str == "spacepoint") return FitterType::SpacepointFitter;
    if(str == "helix") return FitterType::HelixFitter;
    return FitterType::Unknown;
}




Int_t main(Int_t argc, char** argv)
{
    // Options for analysis
    Options opts;
    FitterType fitter = FitterType::Unknown;

    Int_t opt;
    while((opt = getopt(argc, argv, "F:e:M:t:spq")) != -1)
    {
        switch (opt)
        {
            case 'F': 
                fitter = GetFitterTypeFromString(optarg);
                break;
            case 'e':
                opts.processAll = false;
                opts.event = stoi(optarg);
                break;
            case 'M':
                opts.eventMax = stoi(optarg);
                break;
            case 's':
                opts.useSmearing = true;
                break;
            case 't':
                opts.turnMode = true;
                opts.turnID = stof(optarg);
                break;
            case 'p':
                opts.pttrecMode = true;
                break;
            case 'q':
                opts.quietMode = true;
                break;
            case '?':
            default:
            cerr << "\n>>> Usage: " << argv[0]
                 << " -F [planar|spacepoint] -e [eventID] -M [eventMax] -t [nTurns] -s -p -q\n\n"
                 << "\tF: fitter type\n"
                 << "\te: single event mode\n"
                 << "\tM: multiple events mode\n"
                 << "\tt: apply turn analysis\n"
                 << "\ts: apply smearing analysis\n"
                 << "\tp: apply pattern recognition analysis\n"
                 << "\tq: quiet mode\n\n";
            return EXIT_FAILURE;
            
        }
    }

    // Standard error file
    freopen("log.err", "w", stderr);

    // Select the fitter
    switch(fitter)
    {
        case FitterType::PlanarFitter:
            FITALG::PlanarFitter(opts);
            break;
        case FitterType::SpacepointFitter:
            FITALG::SpacepointFitter(opts);
            break;
        case FitterType::HelixFitter:
            FITALG::HelixFitter(opts);
            break;
        default:
            cerr << ">>> Unknown fitter type!" << endl;
            return EXIT_FAILURE;
    }

    // Finally
    return EXIT_SUCCESS;
}