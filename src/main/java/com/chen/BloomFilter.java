package com.chen;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import com.google.common.hash.Hashing;
import java.nio.charset.StandardCharsets;
public class BloomFilter {
    private int n; // size
    private float p; // 误报率
    private int k;
    private BitArray m_bits;

    public BloomFilter(int sz, float FPR, int _k) {
        this.n = sz;
        this.p = FPR;//应该是没用的
        this.k = _k;
        this.m_bits = new BitArray(sz);
    }

    public BitArray getM_bits(){
        return m_bits;
    }

    public void insert(List<Integer> a) {
        for (int idx : a) {
            m_bits.set(idx, true);
        }
    }

    public boolean test(List<Integer> a) {
        for (int idx : a) {
            if (!m_bits.get(idx)) {
                return false;
            }
        }
        return true;
    }

    public void serializeBF(Path BF_file) throws Exception {
        m_bits.serializeBitAr(BF_file);
    }

    public void deserializeBF(Path BF_file) throws Exception {
        m_bits.deserializeBitAr(BF_file);
    }

    public static List<Integer> myhash(String key, int len, int k, int range, int seed) {
        List<Integer> hashvals = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            // Calculate MurmurHash3 hash value
            int hash = Hashing.murmur3_32(seed * k + i)
                    .hashString(key, StandardCharsets.UTF_8)
                    .asInt() % range;
            hashvals.add(hash);
        }
        return hashvals;
    }
}
