
# read references into dictionary
f = open('ictv.list')
text = f.read()
lines = text.splitlines()
names2ids = {name: id for (id,name) in enumerate(lines[1:]) }
f.close()


input = open('ictv.ani')
output = open('ictv.num', 'w')
text = input.read()
lines = text.splitlines()
output.write(lines[0] + '\n')
for line in lines[1:]:
    tokens = line.split('\t')
    ids = [names2ids[tokens[0]], names2ids[tokens[1]]]
    output.write(str(ids[0]) + '\t' + str(ids[1]) + '\t' + tokens[2] + '\n')
input.close()
output.close()
