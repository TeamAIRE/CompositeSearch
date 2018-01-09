/*

        Written by Jananan PATHMANATHAN, 2014-2018
        
        This file is part of CompositeSearch.
        
        CompositeSearch is shared under Creative commons licence: 
        
        Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
        
        See https://creativecommons.org/licenses/by-nc-sa/4.0/

*/

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <map>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <iostream>

// Genes
class geneInfo
{
	public:
		geneInfo();
		geneInfo(unsigned int length, bool composite, unsigned int NbDomains, float noOverlap, unsigned long long int familyID);
		// Set or insert values
		void setLength(unsigned int);
		void setComposite(bool);
		void setNbDomains(unsigned int);
		void setNoOverlap(float);
		void setFamilyID(unsigned long long int);
		void insertComponentFamilyID(unsigned long long int);
		void insertNeighbor(unsigned long long int);
		void sortNeighbors();
		// Get values
		unsigned int getLength();
		bool isComposite();
		unsigned int getNbDomains();
		float getNoOverlap();
		unsigned long long int getFamilyID();
		std::vector<unsigned long long int> getComponentFamiliesID();
		std::vector<unsigned long long int> getNeighbors();
		unsigned long long int getNeighborsSize();
	private:
		unsigned int m_length;
		bool m_composite;
		unsigned int m_NbDomains;
		float m_noOverlap;
		unsigned long long int m_familyID;
		std::vector<unsigned long long int> m_componentFamiliesID;
		std::vector<unsigned long long int> m_neighbors;
};

// EDGES	

class edgeValues
{
	public:
		edgeValues();
		edgeValues(unsigned short int qstart, unsigned short int qend, unsigned short int sstart, unsigned short int send, float pident, double evalue);
		// Set values for the edge
		void setQstart(unsigned short int);
		void setQend(unsigned short int);
		void setSstart(unsigned short int);
		void setSend(unsigned short int);
		void setPident(float);
		void setEvalue(double);
		// Get values 
		unsigned short int getQstart();
		unsigned short int getQend();
		unsigned short int getSstart();
		unsigned short int getSend();
		float getPident();
		double getEvalue();
	private:
		// Query
		unsigned short int m_qstart;
		unsigned short int m_qend;
		// Subject
		unsigned short int m_sstart;
		unsigned short int m_send;
		// Hits values
		float m_pident;
		double m_evalue;
};

// FAMILIES

class familyInfo
{
	public:
		familyInfo();
		familyInfo(float connectivity, bool composite);
		// Set or add values
		void setConnectivity(float);
		void setComposite(bool);
		void insertGeneID(unsigned long long int);
		void sortGenes();
		// Get values
		float getConnectivity();
		bool isComposite();
		std::vector<unsigned long long int> getGenesID();
		unsigned long long int getSize();
	private:
		float m_connectivity;
		bool m_composite;
		std::vector<unsigned long long int> m_genes;
};


edgeValues getEdgeValues(std::map<std::pair<unsigned long long int,unsigned long long int>, edgeValues> &edges, unsigned long long int a, unsigned long long int b);

edgeValues reverseEdgeValues(edgeValues edge);

bool checkTriplet(edgeValues edge1, edgeValues edge2, unsigned int MAX_OVERLAP);

bool FoundIn(unsigned long long int value, std::vector<unsigned long long int>& List);

bool checkCoverage(unsigned long long int node1, unsigned long long int node2, unsigned short int minCov, std::map<unsigned long long int, geneInfo>& genes, std::map<std::pair<unsigned long long int, unsigned long long int>, edgeValues>& edges);

unsigned int computeMedian(std::vector<unsigned int>& val);

bool checkOverlap(unsigned int MAX_OVERLAP, unsigned int inter1Start, unsigned int inter1End, unsigned int inter2Start, unsigned int inter2End);

bool mergeRegion(unsigned int inter1Start, unsigned int inter1End, unsigned int inter2Start, unsigned int inter2End, unsigned short int minCov, unsigned int MAX_OVERLAP);


#endif // FUNCTIONS_H_INCLUDED
