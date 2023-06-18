

a = open("/Users/ksmpooh/Downloads/rosalind_dna.txt","r")

t = a.readline()
t = t.strip()
out = []
for i in ["A","C","G","T"]:
   out.append(t.count(i))

print(out)
