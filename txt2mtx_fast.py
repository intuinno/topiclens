import operator
import sys
import os

do_stemming = False
input_file = sys.argv[1]
output_file = sys.argv[1]

def truncate_word(word):
	start = 0
	while start < len(word) and word[start].isalpha() == False:
		start += 1
	end = len(word)
	while end > start and word[end-1].isalpha() == False:
		end -= 1
	truncated = word[start:end].lower()
	for letter in truncated:
		if letter.isalpha():
			break
	else:
		return ''
	try:
		truncated.decode('ascii')
	except UnicodeDecodeError:
		return ''
	if truncated.find('http://') == 0:
		return ''
	if do_stemming == True:
		if len(truncated) == 0:
			return ''
#		else:
#			return stemmer.stem(truncated, 0, len(truncated)-1)
	else:
		return truncated

def tokenization(text):
	if len(text) == 0:
		return
	start_pos = 0
	while start_pos < len(text):
		while start_pos < len(text):
			if text[start_pos].isalpha():
				break
			else:
				start_pos += 1
		end_pos = start_pos
		while end_pos < len(text):
			if text[end_pos].isalpha():
				end_pos += 1
			else:
				break
		word = text[start_pos:end_pos].lower()
		if word.find('urgent') == -1:
			yield text[start_pos:end_pos]
		start_pos = end_pos

def read_txt(text):
	bag_words = dict()
	for word in tokenization(text):
		truncated = truncate_word(word)
		if truncated != '':
			try:
				bag_words[truncated] += 1
			except KeyError:
				bag_words[truncated] = 1
	return bag_words

stop_list = set()
f_stop = open('english.stop')
for line in f_stop:
	stop_list.add(line[0:-1])
f_stop.close()

if os.path.isdir(input_file):
	bag_words = dict()
	doc_count = 0
	for fname in os.listdir(input_file):
		f_each = open (os.path.join(input_file,fname))
		text = f_each.read()
		f_each.close();
		bag_words_one = read_txt(text)
		for word in bag_words_one:
			try:
				bag_words[word] += bag_words_one[word]
			except KeyError:
				bag_words[word] = bag_words_one[word]
		doc_count += 1
else:
	f_tweets = open(input_file)
	bag_words = dict()
	doc_count = 0
	for line in f_tweets:
		text = line[0:-1]
		bag_words_one = read_txt(text)
		for word in bag_words_one:
			try:
				bag_words[word] += bag_words_one[word]
			except KeyError:
				bag_words[word] = bag_words_one[word]
		doc_count += 1
		print doc_count
	f_tweets.close()

voc_file = open('tmp1/vocabulary.txt', 'w')
word_map = dict()
count = 0
for word in sorted(bag_words):
	if word in stop_list:
		continue
	voc_file.write(word + '\t' + str(bag_words[word]) + '\n')
	word_map[word] = count
	count += 1
voc_file.close()

f_mtx = open('tmp1/tmp.mtx', 'w')
#f_mtx.write('%%MatrixMarket matrix coordinate real general\n%\n')
f_mtx.write('%%MatrixMarket matrix coordinate real general\n')
#f_mtx.write(str(len(word_map)) + ' ' + str(doc_count) + ' 0\n')
doc_count_tot = doc_count
word_map_tot_len = len(word_map)
doc_count = 0
line_count = 0
if os.path.isdir(input_file):
	for fname in os.listdir(input_file):
		f_each = open (os.path.join(input_file,fname))
		text = f_each.read()
		f_each.close();
		bag_words_one = read_txt(text)
		for word in sorted(bag_words_one):
			try:
				word_idx = word_map[word]
				tmp = (str(word_idx+1) + ' ' + str(doc_count+1) + ' ' + str(bag_words_one[word]) + '\n')
				line_count +=1
			except KeyError:
				continue
	print str(word_map_tot_len) + ' ' + str(doc_count_tot) + ' ' + str(line_count) + '\n'
	f_mtx.write(str(word_map_tot_len) + ' ' + str(doc_count_tot) + ' ' + str(line_count) + '\n')
	f_flist = open('tmp1/flist.txt', 'w');
	for fname in os.listdir(input_file):
#		print fname
		f_flist.write(fname + '\n')
		f_each = open (os.path.join(input_file,fname))
		text = f_each.read()
		f_each.close();
		bag_words_one = read_txt(text)
		for word in sorted(bag_words_one):
			try:
				word_idx = word_map[word]
				f_mtx.write(str(word_idx+1) + ' ' + str(doc_count+1) + ' ' + str(bag_words_one[word]) + '\n')
			except KeyError:
				continue
		doc_count += 1
	f_flist.close()
else:
	f_tweets = open(input_file)
	for line in f_tweets:
		text = line[0:-1]
		bag_words_one = read_txt(text)
		for word in sorted(bag_words_one):
			try:
				word_idx = word_map[word]
				f_mtx.write(str(word_idx+1) + ' ' + str(doc_count+1) + ' ' + str(bag_words_one[word]) + '\n')
				line_count +=1
			except KeyError:
				continue
		doc_count += 1
	f_tweets.close()
f_mtx.close()
