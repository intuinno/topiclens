# -*- coding: utf-8 -*-

from flask import Flask, request, g, render_template
from scipy import sparse, io
import numpy as np
from matlab import engine
import os, json


# Configuration
app = Flask(__name__)
app.config['DEBUG'] = True
app.config.from_object(__name__)
app.config.from_envvar('FLASKR_SETTINGS', silent=True)


# Routing
@app.before_first_request
def before__first_request():
	global eng
	global mappedX
	global cl_idx
	global Wtopk
	global voca
	eng = engine.start_matlab()
	eng.cd(os.getcwd()+'\\..\\')
	[mappedX, cl_idx, Wtopk_idx,voca] = eng.main_topic_tsne(nargout=4)

	Wtopk = []
	for idxArray in Wtopk_idx:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk.append(tempArray)

	cl_idx = cl_idx[0]

@app.teardown_request
def teardown_request(exception):
	print('Teardown arose!'.format(exception))


@app.route('/get_subTopic')
def get_subTopic():
	global eng
	global voca
	idx = json.loads(request.args.get('idx'))

	[mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = eng.sub_topic_tsne(idx,nargout=3)

	Wtopk_sub = []
	for idxArray in Wtopk_idx_sub:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk_sub.append(tempArray)

	cl_idx_sub = cl_idx_sub[0]

	return mappedX_sub,cl_idx_sub, Wtopk_sub

# keyword 입력받음
@app.route('/')
def form():
	global mappedX
	global cl_idx
	global Wtopk
	print "rendering..."
	return render_template('tsne.html', mappedX=mappedX, cl_idx=cl_idx, Wtopk= Wtopk)


# Execute the main program
if __name__ == '__main__':
	app.run(host='0.0.0.0',port=5003)
