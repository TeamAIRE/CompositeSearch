/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef COMPOSITEFAMILIES_H_INCLUDED
#define COMPOSITEFAMILIES_H_INCLUDED

#include "functions.h"
#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

void compositeFamilies( std::map<unsigned int, std::list<unsigned long long int> >& familiesToCheck, std::map<unsigned long long int, geneInfo>& genes, std::map<unsigned long long int, familyInfo>& geneFamilies, unsigned int l, std::string timeInfo);

void runCompositeFamilies( std::map<unsigned long long int, geneInfo>& genes, std::map<unsigned long long int, familyInfo>& geneFamilies, unsigned int nCpu, std::string timeInfo);

#endif // COMPOSITEFAMILIES_H_INCLUDED

