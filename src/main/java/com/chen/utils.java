package com.chen;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.io.InputStreamReader;

public class utils {

    public static List<String> getKmers(String filepath) throws IOException {
        String ext = getExtension(filepath);//获取扩展名
        List<String> keys = new ArrayList<>();
        if (ext.equals(".ctx")) {
            keys = getCtxData(Path.of(filepath));
        } else if (ext.equals(".txt") || ext.equals(".out")) {
            keys = getTxtData(filepath);
        } else if (ext.equals(".fna") || ext.equals(".fasta") || ext.equals(".fa") || ext.equals(".ffn")) {
            keys = getTxtData(filepath);
        } else {
//            System.out.println("File extension " + ext + " not recognized! Defaulting to .txt format");
            keys = getTxtData(filepath);
        }
        return keys;
    }

    private static List<String> getTxtData(String txtfile) throws IOException {
        List<String> keys = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(txtfile));
        String line;
        while ((line = br.readLine()) != null) {
            keys.add(line);
        }
        br.close();
        return keys;
    }

//    private static List<String> getCtxData(String ctxfile) throws IOException {
//        List<String> keys = new ArrayList<>();
//        // Execute command to get data from ctxfile
//        // keys = getTxtData(tempDirectoryPath + ctxfile.getFileName().toString());
//        return keys;
//    }

    public static List<String> getCtxData(Path ctxFile) throws IOException {
        List<String> keys = new ArrayList<>();
        // Run the command and read the output
        Process process = Runtime.getRuntime().exec("cortexpy view graph " + ctxFile.toString());
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                // Assuming each line in the output represents a k-mer
                keys.add(line);
            }
        }
        return keys;
    }

    public static String getExtension(String filepath) {
        Path path = Paths.get(filepath);
        String filename = path.getFileName().toString();
        int lastIndex = filename.lastIndexOf('.');
        if (lastIndex != -1) {
            return filename.substring(lastIndex + 1);
        }
        return ""; // Return empty string if extension not found
    }


    public static List<String> readInputFiles(String filePath) throws IOException {
        List<String> inputFiles = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(filePath))) {
            String line;
            while ((line = reader.readLine()) != null) {
                inputFiles.add(line);
            }
        }

        return inputFiles;
    }
}
