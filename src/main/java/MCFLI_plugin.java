import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.*;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

import ij.IJ;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Analyzer;
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.gui.GenericDialog;
import ij.gui.ImageWindow;
import ij.gui.Plot;
import loci.formats.FormatException;
import loci.formats.gui.*;
import trainableSegmentation.*;


public class MCFLI_plugin implements PlugIn {
        //ver 0.6.0  - Major bugfix to Dimitris Plot. Filtering of Dimitris plot by relative error now default.
        //ver 0.5.0b - 128x128 px res for testing 
        //           - got rid of superfluous regression fitting 
        //           - DPlot error now standard error not Standard Dev 
        //           - added sample number and relative standard error diagnostic plots to accompany DPlot
	//ver 0.5.0 Option to select which channel to segment added
        //ver 0.4.0 include choice of bit depth and handling different max intensities during analysis (FYI pixelValueAnalyser doesn't care). 
        //          Uses bit depth to ignore saturated pixels.
        //ver 0.3.0 various bugfixes
        //
        //ver 0.2.0
	//change log inherited from AICA:-all analysis and output can be selected in the dialog box
	//				-lower ch1 and ch2 thresholds introduced for user request (of limited use but hey, wasn't hard to add and will be just as easy to remove in future) 
	//				-Dimitris' plot (DPlot) added as alternative data presentation to 2D histogram
	//				-oneDimensionalAnalysis added. Could replace oneDimensionalCount in future to streamline code.
	//				-lines calling segmentationThresholdAnalysis removed, but method left in code for now just in case it ever becomes useful
	
        
        @Override
	public void run(String arg) {
		
            int numberChannels = 3; //3 channel imaging only for now, user specification is possible later
            int resXY = 128; //default image resolution hard wired to be 1024x1024 for now 
	    
            //Initial setup: User specifies channel that contains independent variable, whether Independent Variable channel requires weka or simple threshold segmentation
            //and asks for threshold values to use.
            GenericDialog setUp = new GenericDialog("Analysis set-up");
	
            //drop down menu for selecting bit depth
            String bitDs [] = {"8", "12", "16"};
            setUp.addChoice("Select Bit Depth", bitDs, bitDs[0]);
            
            //drop down menu for selecting channel roles
            String channelRoles [] = {"Independent variable", "Dependent variable", "Controlled variable"}; 
            String items [] = {"Ch1", "Ch2", "Ch3"};
            
            for (int i=0; i<3; i++) {
                setUp.addChoice(channelRoles[i],items,"Ch" + Integer.toString(i+1));
            }
            
            String segchlchoice [] = {"None", "Ch1", "Ch2", "Ch3"};
            setUp.addChoice("Choose channel for segmentation", segchlchoice, segchlchoice[1]);
            
            
	    setUp.addCheckbox("Weka Segmentation", false);
	    setUp.addNumericField("Weka segmentation probability threshold", 0, 2);
	    setUp.addNumericField("Set lower intensity threshold", 0, 2);
            
            //checkbox group 1
	    int cbRows2 = 2;
	    int cbColumns2 =1;
	    String group2Labels [] = {"Analyse individual channels", "Create Dimitris plots"};
	    String group2Headers [] ={"Analysis and output options"};
	    boolean group2Logic [] = {false, false};
	    setUp.addCheckboxGroup(cbRows2, cbColumns2, group2Labels, group2Logic, group2Headers);
	    setUp.addNumericField("Set relative error threshold (Set a value of 1.0 to include all points)", 0.1, 2);
	    setUp.setOKLabel("Let's do it!");
	    setUp.setCancelLabel("Cancel plugin");
	    setUp.showDialog();
	    if (setUp.wasCanceled()) return;
	    
            int bitDepth;
            int bitChoice = setUp.getNextChoiceIndex();
            switch (bitChoice) {
                case 0:
                    bitDepth = 8;
                    break;
                case 1:
                    bitDepth = 12;
                    break;
                default:
                    bitDepth = 16;
                    break;
            }
            int maxIntensity = (int)Math.pow(2, bitDepth) - 1; //TODO should eventually replace this with max value from actual data
            
            int channelChoices[] = new int [3]; //index is channel role (0-2 --> independent, dependent, controlled), value is channel ID index (i.e. begin at 0)
            for (int i=0; i<3; i++) {
                channelChoices[i] = setUp.getNextChoiceIndex();
            }
            IJ.log("Independent channel index: " + channelChoices[0] + " Dependent channel index: " + channelChoices[1] + " Controlled channel index: " + channelChoices[2]);
            
            boolean segment = true;
            double segThresh = 0;
            double indieThreshold = 0;
            
            int segChl = 0; 
            int segChoice = setUp.getNextChoiceIndex();
            switch (segChoice) {
                case 0:
                    segment = false;
                    break;
                case 1:
                    segChl = 0;
                    break;
                case 2:
                    segChl = 1;
                    break;
                default:
                    segChl = channelChoices[0]; //default set to segment by independent variable
            }
            if (segChoice > 0) {
                segment = setUp.getNextBoolean();
                segThresh = setUp.getNextNumber();
                IJ.log ("weka selection: " + segment + " " + segThresh);
                indieThreshold = setUp.getNextNumber();
                IJ.log ("thresh selection: " + indieThreshold);
            } else {
                setUp.getNextBoolean();// skip through weka checkbox if unneeded 
            }
	    
	    boolean aic = setUp.getNextBoolean();
	    boolean dp = setUp.getNextBoolean();
            
            double filterLim = setUp.getNextNumber();
            IJ.log ("Relative error threshold: " + filterLim);
            
	    IJ.log ("Analysis choices");
            IJ.log ("AIC: " + aic + " DP: " + dp);
	    
	    //User selects *.lif file for analysis
	    OpenDialog setImageFile = new OpenDialog("Select image file");
		String imDir = setImageFile.getDirectory();
	    String imName = setImageFile.getFileName();
	    String imageFileID = imDir + imName;
	    
	    IJ.log("File path: " + imageFileID);
	    	    
	    
	    //If segmentation is required, ask user for classifier location
	    String classFileID = new String();
	    
	    if (segment==true) {
	    	OpenDialog setClassFile = new OpenDialog("Select Weka Segmentation Classifier");
	    	String classDir = setClassFile.getDirectory();
	    	String className = setClassFile.getFileName();
	    	classFileID = classDir + className;
	    		
	    	IJ.log("classifier path: " + classFileID);
	    }
	    
	    //set-up results table to output goodies
	    ResultsTable rt = Analyzer.getResultsTable();
	    if (rt == null) {
	    	rt = new ResultsTable();
	    	Analyzer.setResultsTable(rt);
	    }
	    rt.setPrecision(3);
	
	    	
	    double pixelValues[][] = new double [numberChannels+1][];    	
	    	
	    try {
			pixelValues = pixelValueAnalyser(imageFileID, classFileID, segment, numberChannels, segChl, resXY);
		} catch (IOException e) {
			IJ.error("AICAplugin cannot read image file (IOException)", e.getMessage());
		} catch (FormatException e) {
			IJ.error("AICAplugin cannot read image file (FormatException)", e.getMessage());
		}	
	    IJ.showStatus("Image analysis and segmentation complete, performing post-analysis..."); 
	    
	    
	    //pixelValueFilter will give every pixel class 1 probability = 1 if segment==false
	    IJ.log("Filtering pixel values by Class and/or intensity threshold...");    
	    double class1Values[][] = pixelValuesFilter(pixelValues, numberChannels, maxIntensity, segThresh, indieThreshold, channelChoices[0]);
	      
	    
	    if (aic==true){
	    	double [][] ch1Histogram = oneDimensionalCount(class1Values[0], maxIntensity, 1.0);
	    	Plot histPlot1 = new Plot ("Ch1 pixel value distribution", "Pixel value", "Frequency",ch1Histogram[0], ch1Histogram[1]);
	    	histPlot1.setColor(Color.RED);
	    	ImageWindow.setNextLocation(300,0);
	    	histPlot1.draw();
	    	histPlot1.show();
		
	    	double [][] ch2Histogram = oneDimensionalCount(class1Values[1], maxIntensity, 1.0);
	    	Plot histPlot2 = new Plot ("Ch2 pixel value distribution", "Pixel value", "Frequency",ch2Histogram[0], ch2Histogram[1]);
	    	histPlot2.setColor(Color.GREEN);
	    	ImageWindow.setNextLocation(300,330);
	    	histPlot2.draw();
	    	histPlot2.show();
                
                double [][] ch3Histogram = oneDimensionalCount(class1Values[2], maxIntensity, 1.0);
	    	Plot histPlot3 = new Plot ("Ch2 pixel value distribution", "Pixel value", "Frequency",ch3Histogram[0], ch3Histogram[1]);
	    	histPlot3.setColor(Color.BLUE);
	    	ImageWindow.setNextLocation(300,660);
	    	histPlot3.draw();
	    	histPlot3.show();
                
	    	double ch1Stats[] = arrayColumnStats(class1Values, 0);
		double ch2Stats[] = arrayColumnStats(class1Values, 1);
                double ch3Stats[] = arrayColumnStats(class1Values, 2);
		    
		    rt.incrementCounter();
		    rt.addLabel("Channel 1 pixel value stats");
		    rt.addValue("Median", ch1Stats[3]);
		    rt.addValue("Lower quartile", ch1Stats[2]);
		    rt.addValue("Upper quartile", ch1Stats[4]);
		    rt.incrementCounter();
		    rt.addLabel("Channel 2 pixel value stats");
		    rt.addValue("Median", ch2Stats[3]);
		    rt.addValue("Lower quartile", ch2Stats[2]);
		    rt.addValue("Upper quartile", ch2Stats[4]);
                    rt.incrementCounter();
		    rt.addLabel("Channel 3 pixel value stats");
		    rt.addValue("Median", ch3Stats[3]);
		    rt.addValue("Lower quartile", ch3Stats[2]);
		    rt.addValue("Upper quartile", ch3Stats[4]);
	    }
	    
	    if (dp==true){
                
                double controlFilteredValues [][][] = new double [3][2][]; 
                int totals[] = new int [3];
                int limit2 = 128;
                int limit1 = 1;
                
                for (int i=0; i<class1Values[0].length; i++){
                    if (class1Values[channelChoices[2]][i]>limit2){
                        totals[2]++;
                    } else if (class1Values[channelChoices[2]][i]>limit1) {
                        totals[1]++;
                    } else {
                        totals[0]++;
                    }
                }
                
                for (int i=0; i<3; i++) {
                        controlFilteredValues[i][0] = new double [totals[i]];
                        controlFilteredValues[i][1] = new double [totals[i]];
                        IJ.log("Number of filter " + i + " datapoints:" + totals[i]);
                }
                
                int counting[] = new int [3];
                for (int i=0; i<class1Values[0].length; i++) {   
                    if (class1Values[channelChoices[2]][i]>limit2) {
                        controlFilteredValues[2][0][counting[2]] = class1Values[channelChoices[0]][i];
                        controlFilteredValues[2][1][counting[2]] = class1Values[channelChoices[1]][i];
                        counting[2]++;
                    } else if (class1Values[channelChoices[2]][i]>limit1) {
                        controlFilteredValues[1][0][counting[1]] = class1Values[channelChoices[0]][i];
                        controlFilteredValues[1][1][counting[1]] = class1Values[channelChoices[1]][i];
                        counting[1]++;
                    } else {
                        controlFilteredValues[0][0][counting[0]] = class1Values[channelChoices[0]][i];
                        controlFilteredValues[0][1][counting[0]] = class1Values[channelChoices[1]][i];
                        counting[0]++;
                    }
                }
                
                
	    	double [][][] specialAnalysis = new double [3][][]; 
                for (int i=0; i<3; i++){    
                    specialAnalysis[i] = oneDimensionalAnalysis(controlFilteredValues[i], maxIntensity, 1);
                }
                

                int[] filterTotals = new int [3];
                
                for (int i=0; i<specialAnalysis[0][0].length; i++){
                    if (specialAnalysis[0][10][i]<=filterLim && specialAnalysis[0][10][i]>0){
                            filterTotals[0]++;
                    }
                    if (specialAnalysis[1][10][i]<=filterLim && specialAnalysis[1][10][i]>0){
                            filterTotals[1]++;
                    }
                    if (specialAnalysis[2][10][i]<=filterLim && specialAnalysis[2][10][i]>0){
                        filterTotals[2]++;
                    }
                }
                    
                double[][][] filteredSpecAnalysis = new double [3][3][];
                for (int i=0; i<3; i++){
                    filteredSpecAnalysis[i][0]= new double [filterTotals[i]];
                    filteredSpecAnalysis[i][1]= new double [filterTotals[i]];
                    filteredSpecAnalysis[i][2]= new double [filterTotals[i]];
                }
                    
                int[] filterCounts = new int [3];
                for (int i=0; i<specialAnalysis[0][0].length; i++){
                    if (specialAnalysis[0][10][i]<=filterLim && specialAnalysis[0][10][i]>0){
                        filteredSpecAnalysis[0][0][filterCounts[0]] = specialAnalysis[0][0][i];
                        filteredSpecAnalysis[0][1][filterCounts[0]] = specialAnalysis[0][3][i];
                        filteredSpecAnalysis[0][2][filterCounts[0]] = specialAnalysis[0][9][i];
                        filterCounts[0]++;
                    }
                    if (specialAnalysis[1][10][i]<=filterLim && specialAnalysis[1][10][i]>0){
                        filteredSpecAnalysis[1][0][filterCounts[1]] = specialAnalysis[1][0][i];
                        filteredSpecAnalysis[1][1][filterCounts[1]] = specialAnalysis[1][3][i];
                        filteredSpecAnalysis[1][2][filterCounts[1]] = specialAnalysis[1][9][i];
                        filterCounts[1]++;
                    }
                    if (specialAnalysis[2][10][i]<=filterLim && specialAnalysis[2][10][i]>0){
                        filteredSpecAnalysis[2][0][filterCounts[2]] = specialAnalysis[2][0][i];
                        filteredSpecAnalysis[2][1][filterCounts[2]] = specialAnalysis[2][3][i];
                        filteredSpecAnalysis[2][2][filterCounts[2]] = specialAnalysis[2][9][i];
                        filterCounts[2]++;
                    }
                }   
                    
                Plot specPlot1 = new Plot ("Dimitris Plot", "Independent Variable Intensity", "Mean Dependent Variable Intensity", filteredSpecAnalysis[0][0], filteredSpecAnalysis[0][1]);
                specPlot1.addErrorBars(filteredSpecAnalysis[0][2]);
                specPlot1.draw();
                
                specPlot1.setColor(Color.RED);
                specPlot1.addPoints(filteredSpecAnalysis[1][0], filteredSpecAnalysis[1][1], 1);
                specPlot1.addErrorBars(filteredSpecAnalysis[1][2]);
                
                specPlot1.setColor(Color.BLUE);
                specPlot1.addPoints(filteredSpecAnalysis[2][0], filteredSpecAnalysis[2][1], 1);
                specPlot1.addErrorBars(filteredSpecAnalysis[2][2]);
              
                specPlot1.setLimits(0, maxIntensity, 0, maxIntensity);
                specPlot1.show();
                    
                }
                
                /*if (dpf==false){
                    Plot specPlot1 = new Plot ("Dimitris Plot", "Independent Variable Intensity", "Mean Dependent Variable Intensity", specialAnalysis[0][0], specialAnalysis[0][3]);
                    specPlot1.addErrorBars(specialAnalysis[0][9]);
                    specPlot1.draw();
                
                    specPlot1.setColor(Color.RED);
                    specPlot1.addPoints(specialAnalysis[1][0], specialAnalysis[1][3], 1);
                    specPlot1.addErrorBars(specialAnalysis[1][9]);
                
                    specPlot1.setColor(Color.BLUE);
                    specPlot1.addPoints(specialAnalysis[2][0], specialAnalysis[2][3], 1);
                    specPlot1.addErrorBars(specialAnalysis[2][9]);
              
                    specPlot1.setLimits(0, maxIntensity, 0, maxIntensity);
                    specPlot1.show();
                
                    double maxSamples = 1;
                    for (int i=0; i<3; i++){
                        for (int j=0; j<specialAnalysis[i][8].length; j++){
                            if (specialAnalysis[i][8][j] > maxSamples) {
                                maxSamples = specialAnalysis[i][8][j];
                            }
                        }
                    }
                
                    Plot specPlot2 = new Plot ("Sample numbers", "Independent Variable Intensity", "N", specialAnalysis[0][0], specialAnalysis[0][8]);
                    specPlot2.setLimits(0, maxIntensity, 0, maxSamples);
                    //specPlot2.setLogScaleY(); // remember to set y scale min to 0.1 for log graphs!
                    specPlot2.setColor(Color.RED);
                    specPlot2.addPoints(specialAnalysis[1][0], specialAnalysis[1][8], 1);
                    specPlot2.setColor(Color.BLUE);
                    specPlot2.addPoints(specialAnalysis[2][0], specialAnalysis[2][8], 1);
                    specPlot2.show();
                
                    Plot specPlot3 = new Plot ("Relative standard errors", "Independent Variable Intensity", "Relative standard Error", specialAnalysis[0][0], specialAnalysis[0][10]);
                    specPlot3.setLimits(0, maxIntensity, 0, 1);
                    //specPlot3.setLogScaleY(); // remember to set y scale min to 0.1 for log graphs!
                    specPlot3.setColor(Color.RED);
                    specPlot3.addPoints(specialAnalysis[1][0], specialAnalysis[1][10], 1);
                    specPlot3.setColor(Color.BLUE);
                    specPlot3.addPoints(specialAnalysis[2][0], specialAnalysis[2][10], 1);
                    specPlot3.show();
                    
                } 
            }*/
            
	    IJ.log("I'm finished :)");
	    IJ.showStatus("Analysis finished");
	}
	
	
	public double[][] pixelValueAnalyser(String fileArg, String classArg, Boolean s, int noChls, int sChl, int imSize) throws IOException, FormatException {
		//ver. 1.01 
		//change log:	-correction to segmentation handling so that only one channel of active series is segmented, previously both were processed unnecessarily taking twice the time.
		//				Now new ImagePlus is created from the ImageProcessor of the active series, 
		//			 	-Also incorrect probResult slice selection has been corrected (ImagePlus objects have channel, z-slice and t-frame numbering with base 1)
		
		//open image file and get number of series
		
	    BufferedImageReader imageFile = new BufferedImageReader();
	    imageFile.setId(fileArg);
	    int seriesCount = imageFile.getSeriesCount();
	    IJ.log("File contains " + seriesCount + " series.");
	    double [][] pxVals  = new double [noChls+1][imSize*imSize*seriesCount];
	    
            
	    for (int i=0; i<seriesCount; i++) {
	    	
	    	imageFile.setSeries(i);
	    	int sizeC = imageFile.getSizeC();
	    	int sizeZ = imageFile.getSizeZ();
	    	int sizeX = imageFile.getSizeX();
	    	int sizeY = imageFile.getSizeY();
	    	
	    	IJ.log("C: " + sizeC + "   Z: " + sizeZ + "   X: " + sizeX + "   Y: " + sizeY);
	    	
	    	if (sizeC==noChls && sizeZ==1 && sizeX==imSize && sizeY==imSize) {
	    				
	    		for (int j=0; j<noChls; j++) {
	    		IJ.log("Analysing Series " + (i+1) + "/" + seriesCount + ", Channel " + (j+1));
	    			BufferedImage activeSeries = imageFile.openImage(j);
	    			ImagePlus activeSeriesImP = new ImagePlus();
	    			activeSeriesImP.setImage(activeSeries);	    				    			
	    			
	    			ImageProcessor activeSeriesIP = activeSeriesImP.getProcessor();
                                
	    			float [][] pxValues = activeSeriesIP.getFloatArray();//far faster to use this than getPixels
	    			
	    			for (int k=0; k<sizeY; k++) {
	    				for(int m=0; m<sizeX; m++) {
	    					pxVals [j][(i*sizeX*sizeY)+(k*sizeY)+m] = pxValues[m][k];	    						    					
	    				}	
	    			}
	    			
	    			if (s==true && sChl==j) {
	    				IJ.log("Applying classifier to Channel " + (j+1) + "...");
	    				ImagePlus toBeClassified = new ImagePlus("", activeSeriesIP);
	    				WekaSegmentation segmentator = new WekaSegmentation(toBeClassified);
	    				segmentator.loadClassifier(classArg);
	    				
	    				ImagePlus classResult = segmentator.applyClassifier(toBeClassified, 0, true);
	    				classResult.setSlice(1);
	    				ImageProcessor probIP = classResult.getProcessor();
	    				
	    				float[][] probData = probIP.getFloatArray();
	    				
	    				for (int k=0; k<sizeY; k++) {
		    				for(int m=0; m<sizeX; m++) {		  
		    					pxVals [noChls][(i*sizeX*sizeY)+(k*sizeY)+m] = probData[m][k];
		    				}
	    				}
	    			}
	    			
	    			if (s==false) {
	    				for (int k=0; k<sizeY; k++) {
		    				for(int m=0; m<sizeX; m++) {
		    					pxVals [noChls][(i*sizeX*sizeY)+(k*sizeY)+m] = 1;	    						    					
		    				}
	    				}
	    			}	    		
	    		}	
	    	}
	    }    	
	imageFile.close();
	return pxVals;
	}
	
	
	public double[][] pixelValuesFilter (double  dataColumns[][], int n, int maxI, double sThresh, double indieChlThresh, int indieVarIndex){
		//ver 1.0: Offspring of segmentationClassSeperator in AICA_plugin
                //takes a n-column array and creates a new array from the values in the data columns, filtered by a threshold applied to the values of the last column 
		//
		
		int freqClass1 = 0;
		
		//count number of values in each class based on threshold then initialise arrays to hold them
		for (int i=0; i<dataColumns[0].length; i++){
                    boolean unsaturated = false;
                    for (int j=0; j<n; j++){
                        if (dataColumns[j][i] < maxI) {
                            unsaturated = true;
                        }
                    }
			if (dataColumns[n][i]>= sThresh && dataColumns[indieVarIndex][i]>indieChlThresh && unsaturated) {
				freqClass1++;
			} 
		}		
		
		double filteredClass1[][] = new double [n+1][freqClass1];
		
		int count1=0;
		
		for (int i=0; i<dataColumns[0].length; i++){
                    boolean unsaturated = false;
                    for (int j=0; j<n; j++){
                        if (dataColumns[j][i] < maxI) {
                            unsaturated = true;
                        }
                    }
			if (dataColumns[n][i] >= sThresh && dataColumns[0][i]>indieChlThresh && unsaturated) {
                            for(int j=0; j<n; j++) {
                                filteredClass1[j][count1]= dataColumns[j][i];
                            }
                            count1++;
			}
		}
		return filteredClass1;
	}

	public double[][] oneDimensionalCount(double toBeCounted[], double xMax, double binSize) {
		//calculates weighted frequency counts of a single column of data values within the specified bin size and range
			
			int topBin = (int) (xMax/binSize);
			int totalN = 0;
			double oneDimFreq[][] = new double [2][topBin];
			
			for (int i=0; i<topBin; i++){
				oneDimFreq [0][i] = (i+0.5)*binSize;
			}
			
			for (int i=0; i<toBeCounted.length; i++){
				double assessBin = toBeCounted[i]/binSize;		//divide the data value by the bin size...
				int elementId = (int) assessBin;			//...then casting the double as an integer chops off the decimal places and provides the array element
				
				if (elementId < topBin) {				// add to bin if bin isn't outside of selected range	
					oneDimFreq[1][elementId] = oneDimFreq[1][elementId] + 1;
					totalN++;
				}
			}
			//make into probability distribution rather than frequency count
			for (int i=0; i<topBin; i++){
				oneDimFreq[1][i] = oneDimFreq[1][i]/totalN;
			}
			
			return oneDimFreq;
		}
		
	public double[][] oneDimensionalAnalysis(double toBeCounted[][], double xMax, double binSize) {
		//More complicated version of oneDimensional count:
		//initialises an array of descriptive statistics objects, each descriptive stats object is for a x-value bin.
		//Now in addition to calculating weighted frequency counts of a single column of data values within the specified bin size and range, 
		//can get stats e.g. mean, median etc.
			
			//calculate number of x-bins from the specified max x-value and bin size
			int topBin = (int) (xMax/binSize);
			
			//initialise array of descriptive stats objects
			DescriptiveStatistics [] statsArray= new DescriptiveStatistics [topBin]; 
			for (int i=0; i<topBin; i++){
				statsArray[i] = new DescriptiveStatistics();
			}	
			
			for (int i=0; i<toBeCounted[0].length; i++){
				double assessBin = toBeCounted[0][i]/binSize;		//divide the data value by the bin size...
				int elementId = (int) assessBin;			//...then casting the double as an integer chops off the decimal places and provides the array element
				
				if (elementId < topBin) {				// add to bin if bin isn't outside of selected range	
					statsArray[elementId].addValue(toBeCounted[1][i]);	
				}
			}
			
			//Going to collect all possible stats as they might be useful in future
			//0-8: bin centre; max; min; mean; stdev; 25th percentile; 75th percentile; median; n; std error
			double oneDimStats[][] = new double [11][topBin];
			
			
			for (int i=0; i<topBin; i++){
				
				oneDimStats [1][i] = statsArray[i].getMax();
				oneDimStats [2][i] = statsArray[i].getMin();
				oneDimStats [3][i] = statsArray[i].getMean();
				oneDimStats [4][i] = statsArray[i].getStandardDeviation();
				oneDimStats [5][i] = statsArray[i].getPercentile(25);
				oneDimStats [6][i] = statsArray[i].getPercentile(75);
				oneDimStats [7][i] = statsArray[i].getPercentile(50);
				oneDimStats [8][i] = statsArray[i].getN();
                                oneDimStats [9][i] = oneDimStats [4][i]/Math.sqrt(oneDimStats [8][i]);
                                oneDimStats [10][i] = oneDimStats [9][i]/oneDimStats [3][i];
				
				if (binSize==1){
					oneDimStats [0][i] = i+1;
				} else {
					oneDimStats [0][i] = (i+0.5)*binSize;
				}
			}
			
			return oneDimStats;
		}
	
	public double [][] dataNormalisation (double toBeNormalised[][]){
		double normalised [][] = new double [toBeNormalised[0].length][toBeNormalised[0].length];
		//TODO implement normalisation
		return normalised;
	}
	
	public double[] arrayColumnStats(double toBeCounted[]) {
		//class for calculating statistics of the values in array columns (single column version)		
		DescriptiveStatistics colStats = new DescriptiveStatistics();
		for (int i=0; i<toBeCounted.length; i++) {
			colStats.addValue(toBeCounted[i]);
		}
		double stats[] = new double [7];
		stats[0] = colStats.getMean();
		stats[1] = colStats.getStandardDeviation();
		stats[2] = colStats.getPercentile(25);
		stats[3] = colStats.getPercentile(50);
		stats[4] = colStats.getPercentile(75);
		stats[5] = colStats.getMin();
		stats[6] = colStats.getMax();
		
	return stats;
	}
	
	public double[] arrayColumnStats(double toBeCounted[][], int column) {
		//class for calculating statistics of the values in array columns (accepts arrays with multiple columns, requires column for analysis to be specified)		
		DescriptiveStatistics colStats = new DescriptiveStatistics();
		for (int i=0; i<toBeCounted[column].length; i++) {
			colStats.addValue(toBeCounted[column][i]);
		}
		double stats[] = new double [7];
		stats[0] = colStats.getMean();
		stats[1] = colStats.getStandardDeviation();	
		stats[2] = colStats.getPercentile(25);
		stats[3] = colStats.getPercentile(50);
		stats[4] = colStats.getPercentile(75);
		stats[5] = colStats.getMin();
		stats[6] = colStats.getMax();
		
	return stats;
	}
        
        public static void main (final String... args){
            new MCFLI_plugin().run("");
        }
}