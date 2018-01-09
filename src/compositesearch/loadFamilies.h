/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef LOADFAMILIES_H_INCLUDED
#define LOADFAMILIES_H_INCLUDED

#include "functions.h"
#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

void loadFamilies( std::string geneFamiliesFile, std::map<unsigned long long int, geneInfo>& genes, std::map<unsigned long long int, familyInfo>& geneFamilies, std::map<unsigned long long int, bool>& isVisited, std::string outputDir);
#endif // LOADFAMILIES_H_INCLUDED
