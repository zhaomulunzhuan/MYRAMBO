package com.chen;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Path;

public class BitArray {
    private byte[] A;
    private int ar_size;

    public BitArray(int size) {
        this.ar_size = size;//bit数
        this.A = new byte[(ar_size / 8) + 1];
        for (int i = 0; i < A.length; i++) {
            A[i] = 0; // Clear the bit array
        }
    }

    public byte[] getA(){
        return A;
    }

    public void ANDop(BitArray B) {
        if (A.length != B.A.length) {
            throw new IllegalArgumentException("BitArrays must be of equal length");
        }
        for (int i = 0; i < A.length; i++) {
            A[i] &= B.A[i];
        }
    }

    public boolean empty() {
        for (byte b : A) {
            if (b != 0) {
                return false;
            }
        }
        return true;
    }

    public void serializeBitAr(Path BF_file) throws IOException {
        try (DataOutputStream out = new DataOutputStream(new FileOutputStream(BF_file.toFile()))) {
            out.write(A);
        }
    }

    public int getcount() {
        int count = 0;
        for (byte b : A) {
            count += Integer.bitCount(b & 0xFF);
        }
        return count;
    }

    public void deserializeBitAr(Path BF_file) throws IOException {
        try (DataInputStream in = new DataInputStream(new FileInputStream(BF_file.toFile()))) {
            in.readFully(A);
        }
    }

    public void set(int idx, boolean value) {
        int byteIndex = idx / 8;
        int bitIndex = idx % 8;
        if (value) {
            A[byteIndex] |= (1 << bitIndex); // Set the bit at idx to 1
        } else {
            A[byteIndex] &= ~(1 << bitIndex); // Set the bit at idx to 0
        }
    }

    public boolean get(int idx) {
        int byteIndex = idx / 8;
        int bitIndex = idx % 8;
        return ((A[byteIndex] >> bitIndex) & 1) == 1; // Check if the bit at idx is 1
    }
}
