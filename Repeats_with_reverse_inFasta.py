import itertools

def ReverseComplement(seq):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
   reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
   return reverse_complement

count = 0
with open('//frosty/utkinai2/Project1/Repeats_with_reverse_1.txt', "w") as ReverseFile:
   with open('//frosty/utkinai2/Project1/Repeats_with_reverse_2.txt', "w") as ReverseFile2:
        with open('//frosty/utkinai2/Project1/identified_orphans_repeats.txt', "r") as f:
            for line1, line2 in itertools.zip_longest(*[f] * 2):
                count +=2
                if count < 36000:
                    ReverseFile.write('>r' + line1[1:])
                    ReverseFile.write(ReverseComplement(line2[:-1]))
                else:
                    if count >= 36000:
                        ReverseFile2.write('>r' + line1[1:])
                        ReverseFile2.write(ReverseComplement(line2[:-1]) + '\n')
# count = 0
# for line in open('//frosty/utkinai2/Project1/Repeats_with_reverse_1.txt', "r"):
#     count+=1
#     if count > 49584 and count < 49589:
#         print(line)