/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef DFSCC_H_INCLUDED
#define DFSCC_H_INCLUDED

#include "functions.h"
#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

void computeFamilies( std::map<unsigned long long int, geneInfo>& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, std::map<unsigned long long int, familyInfo>& geneFamilies, std::map<unsigned long long int, bool>& isVisited, unsigned short int minCov, std::string outputDir);

void computeComponentFamilies( std::vector<unsigned long long int>& components, std::map<unsigned long long int, std::vector<unsigned long long int> >& componentFamilies, std::map<unsigned long long int, geneInfo>& genes);
#endif // DFSCC_H_INCLUDED
