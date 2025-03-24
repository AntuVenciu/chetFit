#include <iostream>
#include <unistd.h>

#include "fitteralgorithms.hh"
#include "auxiliaryalgorithms.hh"
#include "options.hh"

using namespace std;


enum class FitterType
{
    PlanarFitter,
    SpacepointFitter,
    Unknown
};



FitterType GetFitterTypeFromString(const std::string& str)
{
    if(str == "planar") return FitterType::PlanarFitter;
    if(str == "spacepoint") return FitterType::SpacepointFitter;
    return FitterType::Unknown;
}




Int_t main(Int_t argc, char** argv)
{
    Options opts;
    FitterType fitter = FitterType::Unknown;

    Int_t opt;
    while((opt = getopt(argc, argv, "F:e:M:t:s")) != -1)
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
                opts.saveMode = true;
                break;
            case 't':
                opts.turnMode = true;
                opts.turnID = stof(optarg);
                break;
            case '?':
            default:
                cerr << "\n>>> Usage: " << argv[0] << " -F [planar|spacepoint] -e [eventID] -M [eventMax] -t [turnID] -s\n";
                return EXIT_FAILURE;
        }
    }


    switch(fitter)
    {
        case FitterType::PlanarFitter:
            FITALG::PlanarFitter(opts);
            break;
        case FitterType::SpacepointFitter:
            FITALG::SpacepointFitter(opts);
            break;
        default:
            cerr << ">>> Unknown fitter type!" << endl;
            return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}