/*

        Written by Jananan PATHMANATHAN, 2014-2016
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef COMPOSITES_H_INCLUDED
#define COMPOSITES_H_INCLUDED

#include "functions.h"
#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

void composites( std::map<unsigned int, std::list<unsigned long long int> >& refSubNodes, std::map<unsigned long long int, geneInfo >& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, std::map<unsigned long long int, familyInfo >& geneFamilies, unsigned short int MAX_OVERLAP, unsigned short int minCov, unsigned int compositeFamMinSize, unsigned int componentFamMinSize, unsigned int l);

void runComposites( std::map<unsigned long long int, geneInfo >& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges, std::map<unsigned long long int, familyInfo>& geneFamilies, std::vector<std::pair<unsigned long long int,unsigned long long int>>& nodesDegree, unsigned short int MAX_OVERLAP, unsigned short int minCov, unsigned int compositeFamMinSize, unsigned int componentFamMinSize, unsigned int nCpu);

#endif // COMPOSITES_H_INCLUDED

