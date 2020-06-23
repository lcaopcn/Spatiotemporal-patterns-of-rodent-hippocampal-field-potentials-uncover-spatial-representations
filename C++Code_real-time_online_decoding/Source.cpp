#include<ipps.h>
#include<ippcore.h>
#include<ippvm.h>
#include"mkl.h"
#include"omp.h"

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cfloat>

#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>

#include <map>
#include <ctime>
#include <cmath>
#include <numeric>
#include <algorithm>

double distcorr(std::vector<std::vector<double> > olebuf);
double distcorr(float ** olebuf, int prsize, int spatial);

const int nminbuf = 3;
const int len = 25;
const int kchan = 128;
const int nspatial = 58;

const int np = 16;
const int nshf = 1000;

const double ripthr = 80;
const double pvalthr = 0.05, scorethr = 13.8115;

int main() {
	// prepare testing data from files
	double temp;

	printf("Loading data ...");
	// initialize buffers
	std::vector<double> ripple(len,0);
	std::vector<double> fripple(len, 0);
	std::vector<std::vector<double> > input(len, std::vector<double>(kchan,0));
		
	// load decoder
	std::vector<std::vector<double> > wdb;
	std::string decoderfile= "wbasis.txt";
	std::ifstream  wdbfile("wbasis.txt");
	if (wdbfile.fail()) {
		printf("\n[Error] wbasis files failed.\n");
		return -1;
	}

	std::string sline;
	getline(wdbfile, sline);
	std::stringstream ss;
	ss << sline;

	double value;
	int SpatialBinNum = 0;
	while (ss >> value) ++SpatialBinNum;
	
	wdb.assign(kchan, std::vector<double>(SpatialBinNum, 0));
	wdbfile.clear();
	wdbfile.seekg(0);
	for (int j = 0; j < kchan; j++)
		for (int ip = 0; ip < nspatial; ip++)
			wdbfile >> wdb[j][ip];
	wdbfile.close();

	// load shuffled decoeder
	std::ifstream  ufile("wdbshfu.txt");
	std::ifstream  rfile("wdbshfr.txt");
	if (ufile.fail() || ufile.fail()) {
		printf("\n[Error] open wdbshf files failed.\n");
		return -1;
	}
	const int shuffle = kchan * nshf;
	std::vector<std::vector<double> > wdbshu(shuffle, std::vector<double>(nspatial, 0));
	std::vector<std::vector<double> > wdbshr(shuffle, std::vector<double>(nspatial, 0));
	for (int i = 0; i < shuffle; i++) {
		for (int j = 0; j < nspatial; j++) {
			ufile >> wdbshu[i][j];
			rfile >> wdbshr[i][j];
		}
	}
	ufile.close(); rfile.close();
	printf("Finished.\n");

	printf("Initialization of intel hilbert tranform buffer ... ");
	// Initialization of intel hilbert tranform buffer
	Ipp32f x[len];
	Ipp32fc y[len];
	
	IppStatus status;
	IppsHilbertSpec* pSpec;
	Ipp8u* pBuffer;
	
	int sizeSpec, sizeBuf;
	
	status = ippsHilbertGetSize_32f32fc(len, ippAlgHintNone, &sizeSpec, &sizeBuf);
	pSpec = (IppsHilbertSpec*)ippMalloc(sizeSpec);
	pBuffer = (Ipp8u*)ippMalloc(sizeBuf);
	
	status = ippsHilbertInit_32f32fc(len, ippAlgHintNone, pSpec, pBuffer);
	printf("Finished.\n");

	// open data streams and output streams
	printf("Open data streaming files ... ");
	std::ifstream rawMUA("raw_FPA.txt"); rawMUA.clear(); rawMUA.seekg(0);
	std::ifstream rawRipple("raw_ripple.txt"); rawRipple.clear(); rawRipple.seekg(0);
	std::ifstream filteredRipple("raw_filtered_ripple.txt"); filteredRipple.clear(); filteredRipple.seekg(0);
	if (rawMUA.fail() || rawRipple.fail() || filteredRipple.fail()) {
		printf("\n[Error] open streaming files failed.\n");
		return -1;
	}
	else
		printf("Succeed.\n");
	remove("out_OLE.txt");  std::ofstream outOLE("out_OLE.txt");
	remove("out_Evaluation.txt");  std::ofstream outEvaluation("out_Evaluation.txt");


	// initializing variables
	bool update = false, evalue = false, reset = false, halt = false, devnt = false;
	double pval = 1, score = 0;
	int ibin = 0;
	std::vector<double> MUAamp(kchan, 0);
	std::vector<double> OLEout(nspatial, 0);
	// ripple amplitude smoothing buffer and parameters
	std::vector<double> ripbuf(3, 0);
	std::vector<std::vector<double> > pr, pru, prr;
	double ripamp = 0, weight[4] = { 0.5417,0.2917,0.1250,0.0417 };

	printf("Starting main loop ... \n");
	// main loop
	std::string line; int itbin = 0;
	double s_initial, s_elapsed, bin_initial, bin_elapsed;
	s_initial = dsecnd(); bin_initial = dsecnd();
	while (std::getline(rawMUA, line)) {
		std::stringstream  lineStream(line);
		
		// read data until one 20ms bin data filled
		if (ibin < len) {
			// read >300Hz filtered neuroal activity data
			double value; int jcol = -1;
			while (lineStream >> value)
				input[ibin][++jcol] = value;

			// read raw and ripple band(140-250Hz) filtered CA1 pyr layer channel data
			rawRipple >> ripple[ibin];
			filteredRipple >> fripple[ibin];

			++ibin;
		}

		if (ibin == len) {
			// ##########################################
			// ### DECODING
			// ##########################################
			++itbin; for (int ip = 0; ip < nspatial; ip++) OLEout[ip] = 0;;
			// Hilbert of ripple signal
			for (int i = 0; i < len; i++)
				x[i] = (Ipp32f)fripple[i];
			status = ippsHilbert_32f32fc(x, y, pSpec, pBuffer);
			ippsMagnitude_32fc((Ipp32fc*)y, x, len);
			ripamp = 0;
			for (int j = 0; j < len; j++)
				ripamp += double(x[j]) / double(len);
			// ripple amplitude smoothing
			if (itbin > 4)
				ripamp = weight[0] * ripamp + weight[1] * ripbuf[0] + weight[2] * ripbuf[1] + weight[3] * ripbuf[2];
			ripbuf[2] = ripbuf[1]; ripbuf[1] = ripbuf[0]; ripbuf[0] = ripamp;

			// Hilbert of lfpMUA and Deoding
			for (int j = 0; j < kchan; j++) {
				MUAamp[j] = 0;
				for (int i = 0; i < len; i++)
					x[i] = (Ipp32f)input[i][j];
				status = ippsHilbert_32f32fc(x, y, pSpec, pBuffer);
				ippsMagnitude_32fc((Ipp32fc*)y, x, len);
				for (int i = 0; i < len; i++)
					MUAamp[j] += double(x[i]) / double(len);
				// matrix multiplication for decoding
				for (int ip = 0; ip < nspatial; ip++) {
					OLEout[ip] += (double)MUAamp[j] * (double)wdb[j][ip];
				}
			}

			if (reset) {
				update = false;  evalue = false; score = 0;
				pr.clear(); pru.clear(); prr.clear();
				reset = false;  devnt = false;
			}

			// event triger
			if (pr.size() < nminbuf) {
				if (ripamp > ripthr && !halt) update = true;
				pval = 1.0; score = 0;
			}
			else {
				update = true;
				evalue = true;
			}

			// update data
			if (update) {
				std::vector<double> pru_tmp(nspatial*nshf, 0);
				std::vector<double> prr_tmp(nspatial*nshf, 0);
				omp_set_num_threads(np);
#pragma omp parallel for 
				for (int ish = 0; ish < nshf; ish++) {
					for (int ip = 0; ip < nspatial; ip++) {
						for (int j = 0; j < kchan; j++) {
							pru_tmp[ish * nspatial + ip] += (double)MUAamp[j] * (double)wdbshu[kchan * ish + j][ip];
							prr_tmp[ish * nspatial + ip] += (double)MUAamp[j] * (double)wdbshr[kchan * ish + j][ip];
						}
					}
				}
				pr.push_back(OLEout); pr.shrink_to_fit();
				pru.push_back(pru_tmp); pru.shrink_to_fit();
				prr.push_back(prr_tmp); prr.shrink_to_fit();
			}
								
			// shuffle evaluation
			if (evalue) {
				double dcorg = distcorr(pr);
				double pvalu = 0, pvalr = 0;

				omp_set_num_threads(np);
#pragma omp parallel for reduction(+: pvalu,pvalr)
				for (int ish = 0; ish < nshf; ish++) {
					float ** pru_eval = new float*[pr.size()];
					float ** prr_eval = new float*[pr.size()];

					for (int i = 0; i < pr.size(); i++) {
						pru_eval[i] = new float[nspatial];
						prr_eval[i] = new float[nspatial];
						for (int ip = 0; ip < nspatial; ip++) {
							pru_eval[i][ip] = pru[i][nspatial * ish + ip];
							prr_eval[i][ip] = prr[i][nspatial * ish + ip];
						}
					}
					
					if (distcorr(pru_eval,pr.size(),nspatial) > dcorg) pvalu += 1.0 / nshf;
					if (distcorr(prr_eval,pr.size(),nspatial) > dcorg) pvalr += 1.0 / nshf;

					for (int i = 0; i < pr.size(); i++) {
						delete[] pru_eval[i]; delete[] prr_eval[i];
					}
					delete[] pru_eval; delete[] prr_eval;
				}
				// monte carlo p value
				pval = pvalu > pvalr ? pvalu : pvalr;

				// pval and score
				if (pval < pvalthr) {
					if (pval < DBL_EPSILON)
						score += -log(0.001);
					else
						score += -log(pval);
				}
				else
					score += 0;

				if (score > scorethr) {
					devnt = true;
					halt = true;
					reset = true;
				}
					
			}

			if (ripamp < ripthr && pval > pvalthr) {
				reset = true; halt = false;
			}
			// ##########################################
			// ### DECODING FINISHED
			// ##########################################
			// after process, refresh data buffers
			ripple.assign(len, 0);
			fripple.assign(len, 0);
			input.assign(len, std::vector<double>(kchan, 0));
			
			ibin = 0;
			s_elapsed = (dsecnd() - s_initial);
			
			printf("\rRipamp: %.3f; Evaluation: %d; pval: %.5f; score: %.5f; Event: %.df; Time used: %.8f", ripamp,
				evalue, pval, score, devnt, s_elapsed*1000); fflush(stdout);
			// data output
			for (int ip = 0; ip < nspatial; ip++)
				outOLE << std::fixed << std::setprecision(8) << OLEout[ip] << " ";
			outOLE << std::endl;
			outEvaluation << std::fixed << std::setprecision(8) << double(itbin)/0.02 << " " << ripamp << " " << evalue << " "
				<< pval << " " << score << " " << devnt << " " << s_elapsed * 1000 << std::endl;

			s_initial = dsecnd();
		}

		
	}
	return 0;
}

double distcorr(std::vector<std::vector<double> > olebuf) {
	const int length = olebuf.size();
	const int nspa = olebuf[0].size();
	std::vector<float> timbin(length, 0);
	std::vector<float> imax(length, 0);
	float meandiff = 0;

	for (int i = 0; i < length; i++) {
		timbin[i] = float(i) + 1;
		imax[i] = float(max_element(olebuf[i].begin(), olebuf[i].end()) - olebuf[i].begin()) + 1;
		if (i > 1) meandiff += fabs(imax[i] - imax[i - 1]);
	}

	float allmean1 = 0, allmean2 = 0;
	std::vector<float> mean1(length, 0), mean2(length, 0);
	std::vector<std::vector<float> > distance1(length, std::vector<float>(length, 0)),
		distance2(length, std::vector<float>(length, 0));
	for (int i = 0; i < length; i++) {
		mean1[i] = 0; mean2[i] = 0;

		for (int j = 0; j < length; j++) {
			distance1[i][j] = fabs(timbin[i] - timbin[j]);
			distance2[i][j] = fabs(imax[i] - imax[j]);

		}
		mean1[i] = std::accumulate(distance1[i].begin(), distance1[i].end(), 0) / float(length);
		mean2[i] = std::accumulate(distance2[i].begin(), distance2[i].end(), 0) / float(length);
		allmean1 += mean1[i] / float(length);
		allmean2 += mean2[i] / float(length);
	}

	float dcov = 0, dvarx = 0, dvary = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			dcov += (distance1[i][j] - mean1[i] - mean1[j] + allmean1) * (distance2[i][j] - mean2[i] - mean2[j] + allmean2);
			dvarx += pow((distance1[i][j] - mean1[i] - mean1[j] + allmean1), 2);
			dvary += pow((distance2[i][j] - mean2[i] - mean2[j] + allmean2), 2);
		}
	}
	dcov /= pow(float(length), 2);
	dvarx /= pow(float(length), 2);
	dvary /= pow(float(length), 2);


 	if (dvarx == 0 || dvary == 0)
		return 0;
	else
		return sqrt(dcov / (sqrt(dvarx) * sqrt(dvary)));
}

double distcorr(float ** olebuf, int prsize, int nspatial) {
	const int length = prsize;
	const int nspa = nspatial;
	float * timbin = new float[length];
	float * imax = new float[length];
	float meandiff = 0;

	for (int i = 0; i < length; i++) {
		timbin[i] = float(i) + 1; float maxval = 0;
		for (int ip = 0; ip < nspa; ip++) {
			if (olebuf[i][ip] > maxval) {
				imax[i] = float(ip) + 1;
				maxval = olebuf[i][ip];
			}
		}
		if (i > 1) meandiff += fabs(imax[i] - imax[i - 1]);
	}

	float allmean1 = 0, allmean2 = 0;
	float * mean1 = new float[length];
	float * mean2 = new float[length];
	float ** distance1 = new float*[length];
	float ** distance2 = new float*[length];
	for (int i = 0; i < length; i++) {
		mean1[i] = 0; mean2[i] = 0;

		distance1[i] = new float[length];
		distance2[i] = new float[length];
		for (int j = 0; j < length; j++) {
			distance1[i][j] = fabs(timbin[i] - timbin[j]);
			distance2[i][j] = fabs(imax[i] - imax[j]);
			mean1[i] += distance1[i][j] / float(length);
			mean2[i] += distance2[i][j] / float(length);
		}
		allmean1 += mean1[i] / float(length);
		allmean2 += mean2[i] / float(length);
	}

	float dcov = 0, dvarx = 0, dvary = 0;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			dcov += (distance1[i][j] - mean1[i] - mean1[j] + allmean1) * (distance2[i][j] - mean2[i] - mean2[j] + allmean2);
			dvarx += pow((distance1[i][j] - mean1[i] - mean1[j] + allmean1), 2);
			dvary += pow((distance2[i][j] - mean2[i] - mean2[j] + allmean2), 2);
		}
	}

	delete[] timbin; delete[] imax;
	delete[] mean1; delete[] mean2;
	for (int i = 0; i < length; i++) {
		delete[] distance1[i]; delete[] distance2[i];
	}
		
	delete[] distance1; delete[] distance2;

	dcov /= pow(float(length), 2);
	dvarx /= pow(float(length), 2);
	dvary /= pow(float(length), 2);

	if (dvarx == 0 || dvary == 0)
		return 0;
	else
		return sqrt(dcov / (sqrt(dvarx) * sqrt(dvary)));
}