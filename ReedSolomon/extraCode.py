########################################################################
########################################################################
###################### Extra code used in CW4 ##########################
########################################################################
########################################################################
import numpy as np
import csv

################################################################
############ testing RS_decode on messages from CW3 ############
################################################################


# message1
p, t = 8, 5
R = np.loadtxt("message1.csv", dtype='int', delimiter=",")
C = RS_decode(R, 2 * t, *table)
print(C)
np.savetxt("message1_decoded.csv", C.reshape(1, len(C)), fmt="%d", delimiter=",")

# message2
p, t = 8, 8
R = np.loadtxt("message2.csv", dtype='int', delimiter=",")
C = RS_decode(R, 2 * t, *table)
print(C)
np.savetxt("message2_decoded.csv", C.reshape(1, len(C)), fmt="%d", delimiter=",")

# message3
p, t = 8, 10
R = np.loadtxt("message3.csv", dtype='int', delimiter=",")
C = RS_decode(R, 2 * t, *table)
print(C)
np.savetxt("message3_decoded.csv", C.reshape(1, len(C)), fmt="%d", delimiter=",")


################################################################
####################### decoding recipe ########################
################################################################


p, t = 8, 16
table = gf_multTable(p)

# recipe decoding
recipe_length = len(np.genfromtxt("recipeC.csv", dtype='int', delimiter=",", usecols=0))
recipe_R = [0]*recipe_length
recipe_CC = recipe_R.copy()
recipe_M = recipe_R.copy()
recipe_ordered = recipe_R.copy()
for i in range(recipe_length):
    recipe_R[i] = np.genfromtxt("recipeC.csv", dtype='int', delimiter=",", skip_header=i, max_rows=1)
    recipe_CC[i] = RS_decode(recipe_R[i], 2*t, *table)
    recipe_M[i] = recipe_CC[i][:len(recipe_CC[i]) - 2 * t]

p, t = 15, 15
table = gf_multTable(p)

# recipe order decoding
order_R = np.genfromtxt("recipeOrderC.csv", dtype='int', delimiter=",", max_rows=1)
order_CC = RS_decode(order_R, 2*t, *table)
order_M = order_CC[:len(order_CC) - 2 * t]

# reorder recipe
for i, j in zip(range(recipe_length), order_M):
    recipe_ordered[i] = recipe_M[j]

with open('recipe.csv', 'w', newline='\n') as f:
    writer = csv.writer(f)
    writer.writerows(recipe_M)

with open('recipeOrder.csv', 'w', newline='\n') as f:
    writer = csv.writer(f)
    writer.writerows([order_M])

with open('recipe_ordered.csv', 'w', newline='\n') as f:
    writer = csv.writer(f)
    writer.writerows(recipe_ordered)

recipe_file = open('recipe_text.txt', 'w')
for i in range(0, recipe_length):
    recipe_file.write(''.join(chr(j) for j in recipe_ordered[i]) + '\n')
recipe_file.close()


################################################################
######################## decoding audio ########################
################################################################


# audio decoding
p = 15
table = gf_multTable(p)

audio_R = np.genfromtxt("audio.csv", dtype='int', delimiter=",", max_rows=1)
writeAudio(audio_R[:len(audio_R) - 2767], 'audioCorrupted.wav')

audio_CC = RS_decode(audio_R, 2767, *table)
audio_M = audio_CC[:len(audio_CC) - 2767]

writeAudio(audio_M, 'audio.wav')


################################################################
####################### decoding part D ########################
################################################################


# secret question decoding
p = 9
table = gf_multTable(p)

message_R = np.genfromtxt("partD.csv", dtype='int', delimiter=",", max_rows=1)
message_CC = RS_decode(message_R, 100, *table)
message_M = message_CC[:len(message_CC) - 100]

file = open('partD.txt', 'w')
file.write(''.join(chr(j) for j in message_M) + '\n')
file.close()


################################################################
#################### part D (decoding date) ####################
################################################################


# secret question
npar = 32
date_R = np.genfromtxt("date.csv", dtype='int', delimiter=",", max_rows=1)
# iterate from 17 down to 8
for p in range(17, 7, -1):
    print(p)
    table = gf_multTable(p)
    date_CC = RS_decode(date_R, 32, *table)
    date_R = date_CC[:len(date_CC) - npar]

file = open('date.txt', 'w')
file.write(''.join(chr(j) for j in date_R) + '\n')
file.close()

################################################################
################################################################
################################################################