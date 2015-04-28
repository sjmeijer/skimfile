/*************************************************
Skimfile creator
Sam Meijer
Adapted from code written by Tom Gilliss

*************************************************/

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TEventList.h>
#include <TEntryList.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TApplication.h>  //This class creates the ROOT Application Environment that interfaces to the windowing system eventloop and eventhandlers

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 	// atof, atoi
#include <iomanip>      // std::setprecision
#include <utility>      // pair<Type1,Type2>
//#include <vector>

using namespace std;

// THIS BLOCK USED TO DETERMINE GHOST CHANNELS //
int mapchannel(int ch) // Input a ROOT channel ch. If it's one of those listed, you'll get a proper cheasy. Otherwise, you get cheasy=1000.
{
    int cheasy = 1000;
	
	// switch(ch){
	// 	case 113:
	// 		chEasy = 0;
	// 		break;
	// 	case 115:
	// }
	//
	
    if(ch == 113) cheasy = 0; //LG channels. Want low gain channels in order to catch upper end of Ga68 B- spectrum 
    if(ch == 115) cheasy = 1;
    if(ch == 119) cheasy = 2;
    if(ch == 121) cheasy = 3;
    if(ch == 145) cheasy = 4;
    if(ch == 147) cheasy = 5;
    if(ch == 149) cheasy = 6;
    if(ch == 151) cheasy = 7;
    if(ch == 153) cheasy = 8;

    if(ch == 112) cheasy = 0; //HG channels 
    if(ch == 114) cheasy = 1;
    if(ch == 118) cheasy = 2;
    if(ch == 120) cheasy = 3;
    if(ch == 144) cheasy = 4;
    if(ch == 146) cheasy = 5;
    if(ch == 148) cheasy = 6;
    if(ch == 150) cheasy = 7;
    if(ch == 152) cheasy = 8;

    return cheasy;
}



// THIS BLOCK USED TO TAG LG VS HG CHANNELS //
int mapchannelLG(int ch) // Input a ROOT channel ch. If it's one of those listed, you'll get a proper cheasy. Otherwise, you get cheasy=1000.
{
    int cheasy = 1000;
    if(ch == 113) cheasy = 0; //LG channels. Want low gain channels in order to catch upper end of Ga68 B- spectrum 
    if(ch == 115) cheasy = 1;
    if(ch == 119) cheasy = 2;
    if(ch == 121) cheasy = 3;
    if(ch == 145) cheasy = 4;
    if(ch == 147) cheasy = 5;
    if(ch == 149) cheasy = 6;
    if(ch == 151) cheasy = 7;
    if(ch == 153) cheasy = 8;
    return cheasy;
}

int main(int argc, char* argv[])
{

	TApplication *App = new TApplication("App", 0, NULL);



	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////// Initializing
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	// ROOT FILES //
	int startrun = 45003001;//45002489;
	int endrun = 45003011;//45002495; // THIS EVENTUALLY NEEDS TO BE UPDATED TO THE LIST OF GOOD RUNS
	int runNumber[128] = {45003001, 45003002, 45003003, 45003004, 45003005};

	// CHANNEL CONSIDERATIONS //
	//int nchannels = 9;

	// TIME CONSIDERATIONS //
	// double timeCorrelationWindow = 1*(67.71)*60; // SSTC window = 1*(Ga68 T_1/2) seconds . Five half-lives may be too many according to Jason's work.
	double startrunStartTime, endrunStopTime; // Used to get the duration (totalTime) of the run
	double totalTime = 0; // The duration of the run. Also used to calculate K-shell event rate

	/*// PLOTS //
	TGraph *graphTimestampComparison = new TGraph(); //("graphTimestampComparison", "", 25713, 1416527900, 1416553613);
	TCanvas *canvTimestampComparison = new TCanvas;*/

	// VARIABLES //
	vector<double>* channel = 0;
	vector<double>* energyCal = 0;
	vector<double>* timestamp = 0;
	double run = 0; 
	double startTime = 0;
	double stopTime = 0;

	// FRIEND BRANCHES and VARIABLES //
	vector<double> globalTimestamp; // There is a vector (globalTimestamp), containing globalTimestampValues of each waveform, for each entry
	double globalTimestampValue; // The values filling each vector globalTimestamp
	double initialTimestamp = 0;
	vector<int> kShellTag; // There is a vector (kShellTag), containing a 0 or 1 for each waveform, for each entry. 0 for non-K-shell, 1 for K-shell
	vector<int> ghostChannelTag; // 0 for real, 1 for ghost
	vector<int> lgChannelTag; // 0 for HG, 1 for LG
	vector<int> sstcTag; // 0 for out of window, 1 for in window
	vector<double> timeSinceEC; // Time between an EC and an event in the SSTC window. For making an exponential decay plot

	// OTHER and FLOW CONTROL //
	double clockFreq = 100000000.0;	// 100 MHz

	char infile[200], infilename[200], friendfile[200], friendfilename[200];
	int nEntriesT, nEntriestF; // number of entries in the original and friend trees
	int initTimestampTag; // This ensures only one event gets stored for each run, in initialTimestampRunArray
	int ch, chEasy, chanSize; 
	int totalWaveformCounter = 0;

	// ENERGY ROI //
	double meanROI, fwhmROI;
	meanROI = 11.0216; //*!*!*!*! Can draw and fit data to get better guesses, later. For now, just giving values. *!*!*!*!
	fwhmROI = 4.99634;

	vector<pair<int,double> > kShellVector; // A vector of pairs (ROOT channel,globalTimestampValue). The "> >" space is important.
	int kShellCounter = 0;
	int kShellVectorSize = 0;


	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////// Tagging Method
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	for(int i=startrun;i<endrun+1;i++)
	{
		// LOAD THE ORIGINAL TREE //
		cout << "----------RUN " << i << "-------------------" << endl;
		sprintf(infile,"mjd_run%d",i);
		sprintf(infilename,"/global/project/projectdirs/majorana/data/mjd/surfprot/data/gatified/P3END/%s.root",infile);
		TChain *t = new TChain("mjdTree"); // Having this inside the for loop prevents files from being linked--the chain just restarts for each run
		t->AddFile(infilename); // Associate ROOT file with the current TChain 
		// ->SetBranchStatus("*",0); // Do not process (0) any (*) of the branches
		// t->SetBranchStatus("channel",1);
		// t->SetBranchStatus("timestamp",1);
		// t->SetBranchStatus("run",1);
		// t->SetBranchStatus("startTime",1);
		// t->SetBranchStatus("stopTime",1);
		// t->SetBranchStatus("energyCal",1);
		// 
			
		t->SetBranchAddress("channel",&channel);
		t->SetBranchAddress("timestamp",&timestamp);
		t->SetBranchAddress("run",&run);
		t->SetBranchAddress("startTime",&startTime);
		t->SetBranchAddress("stopTime",&stopTime);
		t->SetBranchAddress("energyCal",&energyCal);

		nentriest=t->GetEntries(); // Get number of entries in t0   
		cout << "Num Entries this Run = " << nentriest << endl;

		// CREATE A NEW TREE WITH BRANCHES, MIRRORING ORIGINAL TREE'S HIERARCHY/STRUCTURE //
		TTree *tf = new TTree("tf","tf"); // Declare new tree (name, title)
		tf->Branch("globalTimestamp",&globalTimestamp); // Declare a new branch for storing globalTimestamps
		tf->Branch("kShellTag",&kShellTag); // Declare a new branch for storing K-shell tag
		tf->Branch("ghostChannelTag",&ghostChannelTag); // Declare a new branch for storing ghost channel tag
		tf->Branch("lgChannelTag",&lgChannelTag); // Declare a new branch for storing LG channel tag
		tf->Branch("sstcTag",&sstcTag); // Declare a new branch for storing SSTC tag
		tf->Branch("timeSinceEC",&timeSinceEC); // Declare a new branch for storing SSTC tag

		/*
		// A CORRUPT TIMESTAMP DETECTION METHOD //
		TGraph *graphDeltaTimestamp = new TGraph();
		TCanvas *canvDeltaTimestamp = new TCanvas;
		vector<double> initialTimestampArray;
		for(int k=0;k<50;k++) // cover the first 50 events
		{
		t->GetEntry(k); // Read all branches of entry and return total number of bytes read.
		chanSize = channel->size(); // Get the size of the "channel array" corresponding to the entry/event specified above
		for(int j=0; j<chanSize; j++)
		{
		= timestamp->at(j);			
		}
		}
		*/


		initTimestampTag = 0;  // This ensures only one event gets stored for each run, in initialTimestampRunArray	
		for(int k=0;k<nentriest;k++) // loop over each entry in chain (here we're using just one tree at a time). Not accounting for corrupt events at run boundaries. Is this were we should skip initial corrupt events?
		{
			t->GetEntry(k); // Read all branches of entry and return total number of bytes read.
			chanSize = channel->size(); // Get the size of the "channel array" corresponding to the entry/event specified above

			// GET INITIALTIMESTAMPS, startrunStarTime, AND endrunStopTime. TAG GHOST CHANNELS //
			ghostChannelTag.resize(chanSize); // Resize the ghostChannelTag vector to fit however many WFs are in this event
			for(int j=0; j<chanSize; j++) ghostChannelTag[j] = 1; // To start, assume all WFs' channels are ghosts and then determine the nonghosts below
			if(channel->at(0)!=0) // if the first entry in "channel" is not zero, enter loop
			{
				for(int j=0; j<chanSize; j++) // loop through each channel (each entry in the channel branch, which has n entries (n channels))
				{
					ch = channel->at(j); // Get the value of the ROOT chan that is the jth entry 
					cheasy = mapchannel(ch); // Get cheasy from feeding the ROOT channel ch to the mapchannel function
					if(cheasy != 1000) // If ch is not the default 1000, but rather is one of the input ROOT channels
					{
						//// START AT at least k=7 TO SKIP CORRUPT EVENTS //// Can make this more sophisticated later (look at delta from neighbors)
						if(k>6 && initTimestampTag==0)
						{
							initialTimestamp = timestamp->at(j); // Save init stamp for use in calculating globalTimestampValue
							initTimestampTag++; // This ensures that, for each run, only one event gets stored in initialTimestampRunArray
						
							//if(i==startrun) startrunStartTime = initialTimestamp; // This avoids the use of startTime which may be corrupt for some runs
						}
					
						if(i==startrun) startrunStartTime = startTime; // This seems like it's at the mercy of corrupt time stamps; perhaps it should be in the "if k>6" area above.
						if(i==endrun) endrunStopTime = stopTime;

						ghostChannelTag[j] = 0; // Change this tag to a zero b/c the channel is not a ghost
					}
				}
			}		
		
			// CREATE FRIEND TREE WITH NEW VARIABLES //
			globalTimestamp.resize(chanSize); // Clear and resize the vector to suit the size of the current entry
			kShellTag.resize(chanSize);   
			lgChannelTag.resize(chanSize);
			sstcTag.resize(chanSize);
			timeSinceEC.resize(chanSize);                         
			for(int j=0; j<chanSize; j++) // loop through each channel (each entry in the channel branch, which has chanSize entries (chanSize channels))
			{
				ch = channel->at(j); // Get the value of the ROOT chan that is the jth entry 
				//if(k%10000 == 0) cout << "Entry k " << k << " Current channel " << ch << endl;	 		

				globalTimestampValue = startTime + ((timestamp->at(j)-initialTimestamp)/clockFreq);
				globalTimestamp[j] = globalTimestampValue; // Save waveform's globalTimestampValue into the entry's globalTimestamp vector	

				// FIND K-SHELL EVENTS AND ADD TO STACK //
				kShellTag[j] = 0; // Tag a non-K-shell waveform with the "false" bit
				if(mapchannelLG(ch)==1000)
				{
					kShellBounds = channelKShellBounds(ch); // Get the KShell ROI bounds for this particular channel
					if(energyCal->at(j) > kShellBounds.first && energyCal->at(j) < kShellBounds.second)
					{	
						kShellCounter++;
						cout << "K-shell event at entry " << k << " chan " << ch << " || lwrBnd " << kShellBounds.first << " uprBnd " << kShellBounds.second << endl;
						kShellTag[j] = 1; // Tag a K-shell waveform with the "true" bit

						// SPACE SAVING ROUTINE (Not currently working) //
						int kShellVectorSizeEnter = int(kShellVector.size());
						if(kShellVectorSizeEnter>1)
						{
							for(int l=0; l<kShellVectorSizeEnter; l++)
							{
								if((globalTimestampValue - kShellVector[l].second) > timeCorrelationWindow)
								{
									cout << "   kShellVectorSize BEFORE shift " << kShellVectorSize << endl;
									cout << "   kShellVector at l " << std::setprecision(12) << kShellVector[l].second << endl;
									cout << "   kShellVector at l+1 " << std::setprecision(12) << kShellVector[l+1].second << endl;
									for(int m=l; m<(kShellVectorSizeEnter-1); m++)
									{
										kShellVector[m]=kShellVector[m+1]; // Shift slots in vector
									}
									kShellVector.resize(kShellVectorSize-1); // Truncate the vector, cutting off the now empty end slot
									cout << "   kShellVector at l after shift " << std::setprecision(12) <<  kShellVector[l].second << endl;
									kShellVectorSize = int(kShellVector.size()); // for cout purposes
									cout << "   kShellVectorSize AFTER shift " << kShellVectorSize << endl;
								}
							}
						}
						kShellVector.push_back(std::make_pair(ch,globalTimestampValue));
						cout << "kShellCounter " << kShellCounter << endl;
				
						//kShellVectorSize = int(kShellVector.size());
						//cout << "kShell globalTimestampValue " << std::setprecision(12) << globalTimestampValue << endl;
						//cout << "kShellCounter " << kShellCounter << " kShellVectorSize " << int(kShellVector.size()) << " energyCal " << energyCal->at(j) << endl;
					}
				}

				// TAG LG VS HG CHANNELS //
				cheasy = mapchannelLG(ch);
				if(cheasy != 1000) lgChannelTag[j] = 1; 
				else lgChannelTag[j] = 0;

				// TAG SSTC EVENTS //
				sstcTag[j]=0;
				kShellVectorSize = int(kShellVector.size());
				timeSinceEC[j] = -1.0; // Using -1 as a bogus value/placeholder
				if(kShellVectorSize > 0)
				{
					for(int m=0; m<kShellVectorSize; m++)
					{
						// Compare against all K-shell events in the stack
						if((globalTimestampValue-kShellVector[m].second)>0 && (globalTimestampValue-kShellVector[m].second)<timeCorrelationWindow && channel->at(j)==kShellVector[m].first)
						{
							sstcTag[j]=1;
							timeSinceEC[j] = globalTimestampValue-kShellVector[m].second;
						}
					}
				}

				// PLOTS //
				//graphTimestampComparison->SetPoint(totalWaveformCounter,totalWaveformCounter,timestamp->at(j));
				//totalWaveformCounter++;

				// TEST //
				//cout << " k, j " << k << ", " << j << "   globalTimestampValue " << std::setprecision(14) << globalTimestampValue << " globalTimestamp[j] " << globa				lTimestamp[j] << endl;i
	
			}

			tf->Fill(); // Fill the new tree with the new branches/vectors-- one vector for each entry

		}




		// SANITY CHECK OUTPUTS //	
		//nentriestf=tf->GetEntries(); 
		//cout << "Total number of entries in friend tree is " << nentriestf << endl;
		//cout << "startTime "<<std::setprecision(12)<<startTime<<" initialTimestamp "<<std::setprecision(12)<<initialTimestamp<<" intitialTimestamp in sec "<<initialTimestamp/clockFreq<<endl; 

		// TEST TREE STRUCTURE AND FILLING OF GLOBALTIMESTAMP VECTOR //
		//for(int k=0;k<nentriestf;k++) // loop over each entry in chain. Not accounting for corrupt events at run boundaries
		//{
		//        tf->GetEntry(k); // Read all branches of entry and return total number of bytes read.
		//        chanSize = channel->size(); // Get the size of the "channel array" corresponding to the entry/event specified above
		//	for(int j=0; j<chanSize; j++) // loop through each channel (each entry in the channel branch, which has chanSize entries (chanSize channels))
		//	{
		//        	//ch = channel->at(j);
		//		cout << " k, j " << k << ", " << j << endl; //" globalTimestamp " << std::setprecision(12) << globalTimestamp[j] << endl;
		//	}
		//}

		// Don't write over files until hierarchy is good again
		// SAVE tf TO TFile, LET LOOP TAKE YOU TO NEXT FILE //
		//sprintf(friendfile,"friend_mjd_run%d",i);
		//sprintf(friendfilename,"%s.root",friendfile); // /global/u2/g/gilliss/dev/Analysis/PC/SSTC_FriendTrees/gatified/ this is the location I would like them in
		//TFile *friendFile = new TFile(friendfilename,"recreate");
		tf->Write(); // Write the tree to the current directory
		//friendFile->Close();
		//delete t, tf;
	
		//t->AddFriend(tf);
		//t->Scan("channel:lgChannelTag:ghostChannelTag:energyCal:kShellTag:globalTimestamp:sstcTag:timeSinceEC","","precision=14 colsize=14");
		tf->Scan("lgChannelTag:ghostChannelTag:kShellTag:globalTimestamp:sstcTag:timeSinceEC","","colsize=15 precision=14");



	}



	// SANITY CHECK OUTPUTS //
	//cout << "Size of K-shell array is " << kShellCounter << " which should equal " << int(kShellVector.size()) << endl; // these are != after kShellVector efficiency routine
	/*
	cout << "The K-shell array: " << endl;
	kShellVectorSize = int(kShellVector.size());
	for(int h = 0;h<kShellVectorSize;h++)
	{
	cout << "Index " << h << " channel " << kShellVector[h].first << " globalTimestampValue " << std::setprecision(12) << kShellVector[h].second << endl;
	}
	*/

	// FIND K-SHELL RATE // 
	//totalTime = (endrunStopTime - startrunStartTime); // difference in unix timestamps gives totalTime in seconds
	//double kShellRate = (kShellCounter/totalTime); // /2; // /2 to account for LG-HG double counting. There are some events that only get registered in one or the other "gain" though
	//cout << "totalTime in hours " << totalTime/3600 << " Approximate kShellRate = " << kShellRate << endl;

	/*// GLOBALTIMESTAMP VS TIMESTAMP GRAPH // the precision is lacking for some reason
	//for first example draw timestamp vs globatimestamp
	canvTimestampComparison->cd(); 
	char name_TimestampComparison[200];
	sprintf(name_TimestampComparison,"timestamp vs globalTimestamp, Run%d-%d",startrun,endrun); // Dont think this is used here
	graphTimestampComparison->GetYaxis()->SetTitle("timestamp (cycles)");
	graphTimestampComparison->GetXaxis()->SetTitle("globalTimestamp (Unix)"); // 83558358 15
	graphTimestampComparison->GetYaxis()->SetRangeUser(15000000,83600000);
	graphTimestampComparison->GetXaxis()->SetRangeUser(0,totalWaveformCounter);//(1416527998,1416553600); //1416500000
	graphTimestampComparison->SetTitle(name_TimestampComparison);
	graphTimestampComparison->SetMarkerStyle(20);
	graphTimestampComparison->SetMarkerColor(1);
	graphTimestampComparison->Draw("AP");*/


	App->Run();

} // End of main()








































