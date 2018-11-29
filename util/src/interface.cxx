#include "interface.h"
#include "RootTuple.h"

RootTuple *eventHandler;

void rootinit(const char *filename, int lfilename)
{
    eventHandler = new RootTuple(strFtoC(filename, lfilename), "events", "process");
    eventHandler->Initialise();
}


void rootwrite()
{
    eventHandler->Write();
}


void rootclose()
{
    eventHandler->Close();
    delete eventHandler;
}


void rootaddparticle(int barcode, double px, double py, double pz, double e)
{
    eventHandler->AddParticle(barcode, px, py, pz, e);
}


void rootaddevent(double weight)
{
    eventHandler->SetWeight(weight);
    eventHandler->AddEvent();
}


void rootadddouble(double *ptr, const char* branchname, int lbranchname)
{
    eventHandler->SetDoubleBranch(strFtoC(branchname, lbranchname), ptr);
}


void rootaddfloat(float *ptr, const char* branchname, int lbranchname)
{
    eventHandler->SetFloatBranch(strFtoC(branchname, lbranchname), ptr);
}


void rootaddint(int *ptr, const char* branchname, int lbranchname)
{
    eventHandler->SetIntBranch(strFtoC(branchname, lbranchname), ptr);
}


void rootaddbool(bool *ptr, const char* branchname, int lbranchname)
{
    eventHandler->SetBoolBranch(strFtoC(branchname, lbranchname), ptr);
}


void rootaddprocessdouble(double *ptr, const char* branchname, int lbranchname)
{
    eventHandler->AddSingleDouble(strFtoC(branchname, lbranchname), ptr);
}


void rootaddprocessfloat(float *ptr, const char* branchname, int lbranchname)
{
    eventHandler->AddSingleFloat(strFtoC(branchname, lbranchname), ptr);
}


void rootaddprocessint(int *ptr, const char* branchname, int lbranchname)
{
    eventHandler->AddSingleInt(strFtoC(branchname, lbranchname), ptr);
}


void rootaddprocessbool(bool *ptr, const char* branchname, int lbranchname)
{
    eventHandler->AddSingleBool(strFtoC(branchname, lbranchname), ptr);
}


std::string strFtoC(const char *str, int len)
{
    // Converts a FORTRAN string to a C++ string
    int tlen = 0;
    char tem;

    // Counts non-blank characters in a string str until a first blank character or the end of the string is met
    while (str[tlen] != ' ' && (tlen < len) && (tem = str[tlen++], tem));

    char *tstr = new char[tlen+1];
    strncpy(tstr, str, tlen);
    tstr[tlen] = '\0';

    return std::string(tstr);
}

