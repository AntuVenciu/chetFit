#ifndef OPTIONS_HH
#define OPTIONS_HH


struct Options
{
    Bool_t processAll = true;
    Bool_t quietMode = false;
    Bool_t turnMode = false;
    Bool_t useSmearing = false;
    Bool_t pttrecMode = false;
    Bool_t usePrefitter = false;
    Int_t event = -1;
    Int_t eventMax = -1;
    Float_t turnID = -1;
};


#endif  // OPTIONS_H