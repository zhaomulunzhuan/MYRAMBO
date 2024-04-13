package com.chen;

import java.io.IOException;
import java.io.SyncFailedException;
import java.nio.file.Path;
import java.util.List;

public class main {


    public static void main(String[] args) throws Exception {
        String filePath = "D:\\Code\\Idea_Codes\\MYRAMBO_FILE\\MYRAMBO_inputfiles.txt";

        List<String> inputFiles = utils.readInputFiles(filePath);

        int n_per_set = 100000000; //cardinality of each set 数据集基数？布隆过滤器大小？
        float FPR = 0.01F;
        int R_all = 10;
        int B_all = 50;
        int kmersize=31;//根据数据集kmer长度简单写死

//        原始构建并序列化
//        long startbuild=System.nanoTime();
//        RAMBO rambo=new RAMBO(n_per_set,FPR,R_all,B_all,inputFiles);
//
//        long endbuild= System.nanoTime();
//        long buildtime=(endbuild-startbuild)/1_000_000_000;
//        System.out.println("构建时间"+buildtime+"秒");
//
//        String serializeFile="D:\\Code\\Idea_Codes\\MYRAMBO_FILE\\serializeFIle";
//        rambo.serializeRAMBO(Path.of(serializeFile));



        //反序列化
        String serializeFile="D:\\Code\\Idea_Codes\\MYRAMBO_FILE\\serializeFIle";

        long startbuild=System.nanoTime();

        RAMBO rambo=new RAMBO(Path.of(serializeFile));

        long endbuild= System.nanoTime();
        long buildtime=(endbuild-startbuild)/1_000_000_000;
        System.out.println("反序列化构建时间"+buildtime+"秒");


        //查询单个kmer
//        List<String> result=rambo.query("TTTGAGCTAATTAGAGTAAATTAATCCAATC");
//        System.out.println("查询结果");
//        for (String sample:result){
//            System.out.println(sample);
//        }

        long startquery=System.nanoTime();

        //查询一个文件
        String queryfile="D:\\Code\\Idea_Codes\\MYRAMBO_FILE\\query.txt";
        rambo.queryFile(queryfile);

        long endquery=System.nanoTime();
        long querytime=(endquery-startquery)/ 1_000_000;
        System.out.println("查询时间"+querytime+"ms");


    }

}
