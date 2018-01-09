/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#include "mergeFiles.h"
#include <cstdlib>
#include <string>

using namespace std;


void mergeFiles(int nCpu, string method, string fileOut, string timeInfo)
{
	unsigned int i;
	string fileToMerge;
	if( method == "composites" )
	{
		// composites
		system(("echo '#geneID\tLength\tFamilyID\tNbCompFamilies\tNbSingleCompFamilies\tNbDomains\tNoOverLapScore' >> " + fileOut + ".compositesinfo").c_str());
		// coomposite families
		system(("echo '#FamilyID\tSize\tNbcomposites\tPercentComposites\tNbCompFamilies\tNbSingleCompFamilies\tConnectivity\tNbDomains\tNoOverLapScore' >> " + fileOut + ".compositefamiliesinfo").c_str());
		for( i=0; i < nCpu; i++)
		{
			fileToMerge = to_string(i) + ".compositesinfo_" + timeInfo;
			system(("cat " + fileToMerge + " >> " + fileOut + ".compositesinfo").c_str());
			system(("rm " + fileToMerge).c_str());
			fileToMerge = to_string(i) + ".composites_" + timeInfo;
			system(("cat " + fileToMerge + " >> " + fileOut + ".composites").c_str());
			system(("rm " + fileToMerge).c_str());
			fileToMerge = to_string(i) + ".compositefamiliesinfo_" + timeInfo;
			system(("cat " + fileToMerge + " >> " + fileOut + ".compositefamiliesinfo").c_str());
			system(("rm " + fileToMerge).c_str());
			fileToMerge = to_string(i) + ".compositefamilies_" + timeInfo;
			system(("cat " + fileToMerge + " >> " + fileOut + ".compositefamilies").c_str());
			system(("rm " + fileToMerge).c_str());
		}
	}
	else
	{
		for( i=0; i < nCpu; i++)
		{
			fileToMerge = to_string(i) + "." + method + "_" + timeInfo;
			system(("cat " + fileToMerge + " >> " + fileOut + "." + method).c_str());
			system(("rm " + fileToMerge).c_str());
		}
	}
}
//end
