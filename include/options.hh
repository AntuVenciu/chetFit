#ifndef OPTIONS_HH
#define OPTIONS_HH


struct Options
{
    Bool_t processAll = true;
    Bool_t saveMode = false;
    Bool_t turnMode = false;
    Int_t event = -1;
    Int_t eventMax = -1;
    Float_t turnID = -1;
};


#endif  // OPTIONS_H