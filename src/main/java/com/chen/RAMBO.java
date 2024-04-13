package com.chen;


import com.google.common.hash.Hashing;
import net.jpountz.xxhash.XXHash64;
import net.jpountz.xxhash.XXHashFactory;

import java.io.*;
import java.nio.charset.StandardCharsets;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.nio.file.Files;
import java.nio.file.Path;

public class RAMBO {
    private int R;
    private int B;
    private int n;
    private float p;
    private int range;
    private int k;
    private int K;
    private float FPR;
    private BloomFilter[] Rambo_array;
    private List<List<Integer>> metaRambo;
    private List<String> idx_to_name;
    private List<List<Long>> idx_and_r_to_b;
    private Map<String,Integer> name_to_idx;

    public RAMBO(int n, float fpr1, int r1, int b1, List<String> inputFiles) throws IOException {
        FPR = fpr1;
        R = r1;
        B = b1;
        K = inputFiles.size();
        k = 3;
        range = n;
        System.out.println(String.format("Creating RAMBO index with R=%d, B=%d, n=%d, k=%d", R, B, range, k));

        Rambo_array = new BloomFilter[B*R];
        metaRambo = new ArrayList<>(); // 创建一个外层 List
        idx_to_name = new ArrayList<>();
        idx_and_r_to_b = new ArrayList<>();
        name_to_idx = new HashMap<>();

        // 初始化内层 List
        for (int i = 0; i < B * R; i++) {
            metaRambo.add(new ArrayList<>()); // 向外层 List 添加内层 List
        }

        System.out.println("Inserting kmers...");
        for(int b=0;b<B;b++){
            for(int r=0;r<R;r++){
                Rambo_array[B*r+b]=new BloomFilter(range,FPR,k);
            }
        }

        for(int i=0;i<inputFiles.size();i++){
            String cur_file=inputFiles.get(i);
            insertion(cur_file);
        }

        // 输出 name_to_idx
        System.out.println("name_to_idx:");
        for (Map.Entry<String, Integer> entry : name_to_idx.entrySet()) {
            System.out.println(entry.getKey() + ": " + entry.getValue());
        }

        // 输出 idx_to_name
        System.out.println("\nidx_to_name:");
        for (String Name : idx_to_name) {
            System.out.println(Name);
        }

        // 输出 idx_and_r_to_b
        System.out.println("\nidx_and_r_to_b:");
        for (List<Long> hashValues : idx_and_r_to_b) {
            System.out.println(hashValues);
        }

        //检测metaRambo
        System.out.println("metaRambo");
        for(List<Integer> list:metaRambo){
            System.out.println(list);
        }


    }


    public void insertion(String inputFile) throws IOException {
        //获取文件名，不包含扩展名
        Path inputPath = Paths.get(inputFile);
        // 获取文件名（包含扩展名）
        String fileNameWithExtension = inputPath.getFileName().toString();
        // 获取不包含扩展名的文件名
        String fileName = fileNameWithExtension.substring(0, fileNameWithExtension.lastIndexOf('.'));

        List<Long> hashValues=new ArrayList<>();
        if (name_to_idx.containsKey(fileName)) {
            System.err.println(fileName + " already in RAMBO index!");
        } else {
            // 否则，将新的文件名和哈希值添加到 RAMBO 索引中
            hashValues = partitionHash(fileName, fileName.length());//hashValues有R个哈希值，依次对应这个数据集在R次repetition存放的布隆过滤器的分区索引b
            name_to_idx.put(fileName, idx_to_name.size());//数据集名到数据集索引的映射
            idx_to_name.add(fileName);//数据集索引到数据集名的映射
            idx_and_r_to_b.add(hashValues);//根据数据集索引idx，和r，得到数据集idx在repetition r存放的分区索引b
            System.out.println("添加数据集索引"+name_to_idx.get(fileName)+":"+fileName);
        }

        List<String> kmers=utils.getKmers(inputFile);//得到输入文件的kmers列表

        for(int r=0;r<R;r++){//计算当前数据集在repetition r存放的分区b（即布隆过滤器）在metaRambo中的索引，将这个数据集的索引加入到metaRambo[b][r](这样表示，但metaRambo是一维）
            int partitionNum= (int) (B*r+hashValues.get(r));
            metaRambo.get(partitionNum).add(name_to_idx.get(fileName));
        }//metaRambo[b][r]存放的是repetition r中的b分区存储的所有数据集索引

        for(String kmer:kmers){//对于数据集中的每个kmer
            for(int r=0;r<R;r++){
                List<Integer> temp=kmerhash(kmer,kmer.length(),k,range,r);//k个哈希值
                BloomFilter bloomFilter = Rambo_array[(int) (B * r + hashValues.get(r))];//找到这个数据集在repetition r存放的分区b，也就是hashValues.get(r)，将kmer存放进去
                bloomFilter.insert(temp);

            }
        }

    }

    public List<String> query(String querykey){//查询一个kmer
        int idx_to_name_size = idx_to_name.size();
        BitArray bitarray_K = new BitArray(idx_to_name_size); // 创建一个位数组，大小为 idx_to_name 的大小，即当前存入的数据集数量
        List<Integer> check = kmerhash(querykey, querykey.length(), k, range, 0); // 查询关键字进行 k 次哈希，check 为 k 个哈希索引
        Set<Integer> toCheckNext = new HashSet<>();
        int bfHitCount = 0; // 命中的布隆过滤器数量
        int potentialSampleCount = 0; // 候选数据集数量
        for(int b=0;b<B;b++){
            if (Rambo_array[b].test(check)){//如果Rambo_array[b]这个布隆过滤器包含查询元素
//                System.out.println(b+"分区包含查询元素");
                for (int idx: metaRambo.get(b)){//遍历当前命中分区存储的数据集索引
                    toCheckNext.add(Math.toIntExact(idx_and_r_to_b.get(idx).get(1)));
//                    System.out.println("repetition 1需要检查的数据集所在分区"+toCheckNext);
                    //idx是repetition 0命中的分区b中的一个数据集的索引，idx_and_r_to_b[idx][1]是这个数据集在repetition 1存放的分区号，
                    //作为下一个repetition要检查的partition之一（这样做可以减少不必要的查询）
                    bitarray_K.set(idx,true);//在bitarray_K中将数据集索引idx对应的位置设置为1
                }
                bfHitCount+=1;
            }
        }
        potentialSampleCount=bitarray_K.getcount();
//        System.out.println( "r=0: 命中的布隆过滤器数="+bfHitCount+"\t 候选数据集数量="+potentialSampleCount+"\t\t\t");
        bfHitCount=0;
        for(int r=1;r<R;r++){//检索repetition 1-（R-1）
            Set<Integer> to_check = new HashSet<>(toCheckNext); // 创建 to_check 的副本
//            System.out.println(to_check); // 打印副本内容
            toCheckNext.clear(); // 只清空 toCheckNext
//            System.out.println(to_check); // 副本内容不受影响
            List<Integer> checkindex=kmerhash(querykey,querykey.length(),k,range,r);
            BitArray bitArrayK1 = new BitArray(idx_to_name.size());
            for (Integer b : to_check) {// 依次检查这个 repetition 下的候选分区
//                System.out.println("检查r="+r+"b="+b);
                if (Rambo_array[b + B * r].test(checkindex)) {// 如果 ramboArray[b + B * r] 这个分区包含查询元素
//                    System.out.println("r="+r+"b="+b+"包含查询元素");
                    for (Integer idx : metaRambo.get(b + B * r)) {// 遍历当前命中分区存储的数据集索引
                        if (r < R - 1) {// 还没检查完 R 个 repetition
                            toCheckNext.add(Math.toIntExact(idx_and_r_to_b.get(idx).get(r + 1)));// toCheckNext 存放下一个 repetition 要查询的分区
                        }
                        bitArrayK1.set(idx, true);// bitArrayK1 中将可能包含查询元素的数据集索引对应的位置设置为 1
                    }
                    bfHitCount++;
                }
            }
            bitarray_K.ANDop(bitArrayK1);// bitArrayK 是 repetition 0 中查询得到的可能包含查询元素的数据集，bitArrayK1 是 repetition 1-(R-1) 的查询结果
            // bitarray_K 是一个长度为数据集数量的位数组，此时 bitArrayK 中为 1 的索引对应的数据集就是查询结果
        }
        // 获取候选数据集数量
        int potential_sample_count = bitarray_K.getcount();
//        System.out.println("最终: 命中的sample_count=" + potential_sample_count);

        // 创建一个字符串列表，用于存储命中的数据集名称
        List<String> ret_samples = new ArrayList<>();

        // 在位数组中寻找值为1的位置，然后将命中的数据集名称添加到结果列表中
        int idx = 0;
        while ((idx = findNextSetBit(bitarray_K,idx)) != -1) {
//            System.out.println(idx);
            ret_samples.add(idx_to_name.get(idx));
            idx++;
        }

        // 返回命中的数据集名称列表
        return ret_samples;

    }

    public void queryFile(String filePath){//一个文件中有多个查询长序列，查询每一个并把查询结果写入输出文件
        String queryresultFile = "D:/Code/Idea_Codes/MYRAMBO_FILE"+"/"+"query_result.txt";//存放查询结果
        try(
                BufferedReader reader=new BufferedReader(new FileReader(filePath));
                BufferedWriter writer=new BufferedWriter(new FileWriter(queryresultFile))){
            String line;
            String sequence="";
            while ((line=reader.readLine())!=null){
                if(line.startsWith(">")){
                    //查询
                    if (!sequence.isEmpty()){
                        writer.write(sequence+"\n");
                        querySequence(sequence,writer);
                        writer.write(line+"\n");
                    }else {
                        writer.write(line+"\n");
                    }
                    sequence="";
                }else {
                    sequence+=line.trim().toUpperCase();
                }
            }
            if(!sequence.isEmpty()){
                writer.write(sequence + "\n");
                //查询最后一段序列
                querySequence(sequence,writer);
            }
        }catch (IOException e){
            System.err.println(e);
        }
    }
    public void querySequence(String sequence, BufferedWriter writer) throws IOException {//查找长序列，每个kmer都存在才报告序列存在
        int kmersize=31;//根据数据集kmer长度简单写死
        List<String> kmerList=new ArrayList<>();
        // 切割sequence并将长度为kmersize的子字符串加入kmerList
        for (int i = 0; i <= sequence.length() - kmersize; i++) {
            String kmer = sequence.substring(i, i + kmersize);
            kmerList.add(kmer);
        }

        List<String> result=new ArrayList<>(query(kmerList.get(0)));
        for(String kmer:kmerList){
            result.retainAll(query(kmer));
        }

        writer.write("查询结果\n");
        // 将查询结果写入到结果文件
        if (!result.isEmpty()){
            for (String datasetName : result) {
                writer.write(datasetName + "\n");
            }
        }else {
            writer.write("未查询到包含查询序列的数据集"+"\n");
        }

    }

    public List<Long> partitionHash(String key, int len) {//分区，进行R次哈希，得到每个数据集在每个repetition存储的分区号b
        List<Long> hashValues = new ArrayList<>();

        // Iterate over k hash functions
        for (int i = 0; i < R; i++) {
            // Calculate hash using MurmurHash algorithm
            long hash = Hashing.murmur3_128(i).hashString(key, StandardCharsets.UTF_8).asLong();

            // Map hash value to the specified range
            long mappedHash = Math.abs(hash) % B;

            hashValues.add(mappedHash);
        }

        return hashValues;
    }

    public static List<Integer> kmerhash(String key, int len, int k, int range, int seed) {//对kmer计算k个哈希值
        List<Integer> hashvals = new ArrayList<>();

        XXHashFactory factory = XXHashFactory.fastestInstance();
        XXHash64 hash64 = factory.hash64();

        for (int i = 0; i < k; i++) {
            long hash = hash64.hash(
                    key.getBytes(StandardCharsets.UTF_8),
                    0,
                    key.length(),
                    seed + i * seed
            );

            // 对哈希值取模得到范围内的整数，并将其添加到结果列表中
            hashvals.add((int) (Math.abs(hash) % range));
        }

        return hashvals;
    }

    // 辅助方法：在位数组中找到下一个值为1的位置
    public int findNextSetBit(BitArray bitarray, int startIdx) {
        int idx = startIdx;
        int arrayLength = bitarray.getA().length;
        while (idx < arrayLength * 8) {
            int byteIndex = idx / 8;
            int bitIndex = idx % 8;
            if (((bitarray.getA()[byteIndex] >> bitIndex) & 1) == 1) {
                return idx;
            }
            idx++;
        }
        return -1; // 没有找到值为1的位置
    }

    public RAMBO(Path rambo_dir) throws Exception {//根据序列化文件创建RAMBO
        String line;
        name_to_idx = new HashMap<>();
        idx_to_name = new ArrayList<>();
        metaRambo = new ArrayList<>();
        idx_and_r_to_b=new ArrayList<>();

        // Read idx_to_name.txt
        BufferedReader idxReader = new BufferedReader(new FileReader(rambo_dir.resolve("idx_to_name.txt").toFile()));
        int idx = 0;
        while ((line = idxReader.readLine()) != null) {
            if (!line.isEmpty()) {
                name_to_idx.put(line, idx);
                idx_to_name.add(line);
                idx_and_r_to_b.add(new ArrayList<>());
                idx++;
            }
        }
        idxReader.close();

        // Read metarambo.txt
        BufferedReader metaReader = new BufferedReader(new FileReader(rambo_dir.resolve("metarambo.txt").toFile()));
        line = metaReader.readLine(); // Read the first line
        String[] parts = line.split(" ");
        R = Integer.parseInt(parts[0]);
        B = Integer.parseInt(parts[1]);
        range = Integer.parseInt(parts[2]);
        k = Integer.parseInt(parts[3]);
        float FPR= 0.01F;//没用

        // Initialize Rambo_array
        Rambo_array = new BloomFilter[B * R];
        for (int i = 0; i < B * R; i++) {
            Rambo_array[i] = new BloomFilter(range,FPR,k);
            metaRambo.add(new ArrayList<>());
        }

        int r = 0, b = 0;
        while ((line = metaReader.readLine()) != null) {
            if (line.startsWith("#")) {
                line = line.substring(2); // Remove the "#" and space
                parts = line.split(" ");
                r = Integer.parseInt(parts[0]);
                b = Integer.parseInt(parts[1]);
            } else {
                List<Integer> metaList = metaRambo.get(B * r + b);
                int sampleIdx = Integer.parseInt(line);
                metaList.add(sampleIdx);
                idx_and_r_to_b.get(sampleIdx).add((long) b);
            }
        }
        metaReader.close();

        deserializeRAMBO(rambo_dir);

        // 输出 name_to_idx
        System.out.println("name_to_idx:");
        for (Map.Entry<String, Integer> entry : name_to_idx.entrySet()) {
            System.out.println(entry.getKey() + ": " + entry.getValue());
        }

        // 输出 idx_to_name
        System.out.println("\nidx_to_name:");
        for (String Name : idx_to_name) {
            System.out.println(Name);
        }

        // 输出 idx_and_r_to_b
        System.out.println("\nidx_and_r_to_b:");
        for (List<Long> hashValues : idx_and_r_to_b) {
            System.out.println(hashValues);
        }

        //检测metaRambo
        System.out.println("metaRambo");
        for(List<Integer> list:metaRambo){
            System.out.println(list);
        }
    }

    public void serializeRAMBO(Path dir) throws Exception {
        System.out.println("Saving RAMBO to " + dir.toString());
        Files.createDirectories(dir);

        // Serialize idx_to_name
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(dir.resolve("idx_to_name.txt").toString()))) {
            for (String name : idx_to_name) {
                writer.write(name + System.lineSeparator());
            }
        }

        try (BufferedWriter metaOut = new BufferedWriter(new FileWriter(dir.resolve("metarambo.txt").toFile()))) {
            metaOut.write(R + " " + B + " " + range + " " + k);
            metaOut.newLine(); // 写入RAMBO的元信息到metarambo.txt

            double totalBits = 0; // 记录R*B个布隆过滤器中置为1的bit总数
            for (int r = 0; r < R; r++) {
                Path repDir = dir.resolve("repitition_" + r);
                Files.createDirectories(repDir); // 每个repetition生成一个目录(文件夹)
                for (int b = 0; b < B; b++) {
                    Path bloomFilterPath = repDir.resolve("filter_" + b + ".bloom"); // 每个partition对应的布隆过滤器生成对应文件
                    totalBits += Rambo_array[b + B * r].getM_bits().getcount();// 累加当前布隆过滤器中置为1的bits
                    Rambo_array[b + B * r].serializeBF(bloomFilterPath); // 将布隆过滤器对应的位数组序列化到对应的bloom_filter_path
                    metaOut.write("# " + r + " " + b);
                    metaOut.newLine(); // 每个r下的每个b写入metarambo.txt
                    for (int sampleIdx : metaRambo.get(B * r + b)) { // 遍历当前r 当前b对应的样本索引
                        metaOut.write(String.valueOf(sampleIdx));
                        metaOut.newLine(); // 写入metarambo.txt
                    }
                }
            }

            System.out.println("Average RAMBO density is " + (totalBits / (B * R * range))); // 计算RAMBO的平均密度
        } catch (IOException e) {
            e.printStackTrace();
        }


        System.out.println("RAMBO serialization complete.");
    }

    public void deserializeRAMBO(Path dir) throws Exception {
        for (int r = 0; r < R; r++) {
            for (int b = 0; b < B; b++) {
                Path repDir = dir.resolve("repitition_" + r);
                Path bloomFilterPath = repDir.resolve("filter_" + b + ".bloom");
                Rambo_array[b + B * r].deserializeBF(bloomFilterPath);
            }
        }

        System.out.println("RAMBO deserialization complete.");
    }


}
