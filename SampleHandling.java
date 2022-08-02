
import java.io.File;
import java.util.*;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;


/**
* preprocesseing the files of the Malicia folder
*/
public class SampleHandling{
	
	HashMap<String, Integer> encodingMap = new HashMap<String, Integer>();
	ArrayList<ArrayList<String>> testSample = new ArrayList<ArrayList<String>>();
	
	/*
	* Read all contents under Malicia/folder 
	* save each sample-file content into an ArrayList
	* save all 1000 sample-files' contents into another ArrayList 
	*/
	public ArrayList<ArrayList<String>> readFiles(String folderPath, int initBound, int sampleSize) throws FileNotFoundException {
		File folder = new File(folderPath);

		File[] files = folder.listFiles();
		System.out.println("There are " + files.length + " malware files in folder " + folderPath);
		// the arrayList saved all opcodes in one sample
		ArrayList <ArrayList<String>> listOfFiles = new ArrayList<ArrayList<String>>();
		
		for (int i = initBound; i < sampleSize; i++) {
			ArrayList<String> fileContent = new ArrayList<String>();

			String filename = files[i].getName();
			int index = filename.lastIndexOf('.');
			
			if (index > 0) {
				String extension = filename.substring(index+1);
				if (extension.equals("txt")) {
					Scanner scan = new Scanner(files[i]);
			    	while (scan.hasNextLine()) {
				    	// save all opcodes in an arrayList
			    		fileContent.add(scan.nextLine());
					}
					listOfFiles.add(fileContent);
					System.out.println("size of file " + i + " - " + files[i].getName() + " is " + fileContent.size());
				}
			}
		}
		System.out.println("Read contents of " + listOfFiles.size() + " files");
		return listOfFiles;
	}

    /**
    * Count occurence of each opcode and save the pair in the HashMap
    */
	public HashMap<String,Integer> countOccurence(ArrayList<ArrayList<String>> contents) {
		HashMap<String,Integer> occurencePair = new HashMap<String, Integer>();
         
		for (int i = 0; i < contents.size(); i++) {
			for (int j = 0; j < contents.get(i).size(); j++) {
				occurencePair.put(contents.get(i).get(j), occurencePair.getOrDefault(contents.get(i).get(j),0) + 1);
			}
		}
		return occurencePair;
	}
    
    /**
	* Sort opcodes and its occurence value based on the frequncey of the occurence from the least to the most
	* assign a number to each opcode according to its occurence
	* the returned hashmap would be the code book used to assign numbers to all the samples
	* reference: https://www.geeksforgeeks.org/sorting-a-hashmap-according-to-values/
	*/
	public HashMap<String, Integer> encodingPair(HashMap<String, Integer> countMap) {
		List<Map.Entry<String, Integer>> list = new ArrayList<Map.Entry<String, Integer>>(countMap.entrySet());
		
		Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {
			public int compare(Map.Entry<String, Integer> obj1, Map.Entry<String, Integer> obj2) {
				return (obj1.getValue()).compareTo(obj2.getValue());
			}
		}); 

		//assign a number to an opcode based on the sorted occurance from the least to the most
		//key: opcode value: an integer
		HashMap<String, Integer> codeBook = new HashMap<String, Integer>();
		for (int i = 0; i < list.size(); i++) {
			codeBook.put(list.get(i).getKey(), i);
		}
		System.out.println("Opcode Encoding: " + codeBook);
		return codeBook;
	}

	/**
	* Generating observation sequences
	*/
	public ArrayList<ArrayList<Integer>> encodingSample(HashMap<String, Integer> map, ArrayList<ArrayList<String>> arrs) {
		ArrayList<ArrayList<Integer>> encodedSampleArr = new ArrayList<ArrayList<Integer>>();

		for (int i = 0; i < arrs.size(); i++) {
			ArrayList<Integer> sample = new ArrayList<Integer>();
			for (int j = 0; j < arrs.get(i).size(); j++) {
				if (map.get(arrs.get(i).get(j)) == null) {
					sample.add(0);
				} else {
					sample.add(map.get(arrs.get(i).get(j)));
				}
			}
			encodedSampleArr.add(sample);
		}
		return encodedSampleArr;
	}

	/**
	* An API to get the observation sequence
	*/
	public ArrayList<ArrayList<Integer>> getObservations(String folderPath, int initBound, int sampleSize) throws FileNotFoundException, IOException {
		ArrayList<ArrayList<String>> fileContent = this.readFiles(folderPath, initBound, sampleSize);
		HashMap<String, Integer> occurencePair = this.countOccurence(fileContent);
		int numOfOpcodes = this.numberOfOpcodes(occurencePair);
		this.encodingMap = this.encodingPair(occurencePair);
		ArrayList<ArrayList<Integer>> obsSeq = this.encodingSample(this.encodingMap, fileContent);
		//this.writeSamples(obsSeq, fileName);
		//System.out.printf("#obs = %d, #opcodes = %d, saved to %s\n", fileContent.size(), numOfOpcodes);
		//System.out.printf("#obs = %d, #opcodes = %d, saved to %s\n", fileContent.size(), numOfOpcodes, fileName);
		return obsSeq;
	}

	/**
	* An API to get test observation sequence
	*/
	public ArrayList<ArrayList<Integer>> getTestObservation(String folderPath, int initBound, int sampleSize) throws FileNotFoundException, IOException {
		ArrayList<ArrayList<String>> fileContent = this.readFiles(folderPath, initBound, sampleSize);
		ArrayList<ArrayList<Integer>> obsSeq = this.encodingSample(this.encodingMap, fileContent);
		System.out.println(fileContent.size());
		System.out.println(obsSeq.size());
		return obsSeq;
	}

	/**
	* An API to get numbers of unique opcode
    * this will be an API to determine M
	*/
	public int getNumUniqueOpcode() {
		return this.encodingMap.size();
	}

	/**
	* Count number of unique opcodes in the countOccurence pair hashmap
	*/
	public int numberOfOpcodes(HashMap<String, Integer> map) {
		return map.size();
	}

	/**
	* Write out the observation sequence file
	* contains opcodes associated numbers
	*/
	public void writeSamples(ArrayList arr, String filePath) throws IOException {
		FileWriter writer = new FileWriter(filePath);
		String lines = "";
		for (int i=0; i< arr.size(); i++) {
			lines = lines.concat(arr.get(i).toString() + "\n");
		}
		writer.write(lines);
		writer.close();
	}

	/*
	* Main Function
	* Used when want to write all encoded opcodes into files
	*/
	public static void main(String[] args) throws FileNotFoundException, IOException {
		SampleHandling obj = new SampleHandling();
		String folderPath1 = "./Malicia (.txt Opcodes)/zbot";
		String folderPath2 = "./Malicia (.txt Opcodes)/winwebsec";
	    String folderPath3 = "./Malicia (.txt Opcodes)/zeroaccess";

		// file contents
		ArrayList<ArrayList<String>> zbotFileContent = obj.readFiles(folderPath1,0, 10);
		// ArrayList<String> winFileContent = obj.readFiles(folderPath2);
		// ArrayList<String> zeroFileContent = obj.readFiles(folderPath3);

		// opcodes and occurence pair
		HashMap<String, Integer> zbotOccurencePair = obj.countOccurence(zbotFileContent);
		// HashMap<String, Integer> winOccurencePair = obj.countOccurence(winFileContent);
		// HashMap<String, Integer> zeroOccurencePair = obj.countOccurence(zeroFileContent);
		
		// numbers of uinque opcodes
		int numOfZbotOpcodes = obj.numberOfOpcodes(zbotOccurencePair);
		// int numOfWinOpcodes = obj.numberOfOpcodes(winOccurencePair);
		// int numOfZeroOpcodes = obj.numberOfOpcodes(zeroOccurencePair);

		// assign a number to an opcode
		HashMap<String, Integer> assignedZbot = obj.encodingPair(zbotOccurencePair);
		// HashMap<String, Integer> assignedWin = obj.encodingPair(winOccurencePair);
		// HashMap<String, Integer> assignedZero = obj.encodingPair(zeroOccurencePair);

		// final encoded observation sample
		ArrayList<ArrayList<Integer>> zbotSample = obj.encodingSample(assignedZbot, zbotFileContent);
		// ArrayList winSample = obj.encodingSample(assignedWin, winFileContent);
		// ArrayList zeroSample = obj.encodingSample(assignedZero, zeroFileContent);

		// write into an observation sequence file
		//obj.writeSamples(zbotSample, "./zbotSamples/");
		// obj.writeSamples(winSample, "./winwebsecSequence.txt");
		// obj.writeSamples(zeroSample, "./zeroSequence.txt");
	}
}