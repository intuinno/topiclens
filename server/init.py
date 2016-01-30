# -*- coding: utf-8 -*-

from werkzeug.contrib.profiler import ProfilerMiddleware
from flask import Flask, request, g, render_template
from flask.ext.triangle import Triangle
from scipy import sparse, io
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
from matlab import engine
import os, json, time
from flask.ext.cors import CORS
from jinja2 import Environment


# Configuration
app = Flask(__name__, static_path='/static')	
Triangle(app)
CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

app.config['DEBUG'] = True
app.config.from_object(__name__)
app.config.from_envvar('FLASKR_SETTINGS', silent=True)


app.config['PROFILE'] = True
app.wsgi_app = ProfilerMiddleware(app.wsgi_app, restrictions=[10])
# @app.after_request
# def after_request(response):
#   response.headers.add('Access-Control-Allow-Origin', '*')
#   response.headers.add('Access-Control-Allow-Headers', 'Content-Type,Authorization')
#   response.headers.add('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE')
#   return response



def supervisedTSNE(distanceMatrix, cl_idx, sameTopicWeight=0.9, differentTopicWeight=1.1):
	n = distanceMatrix.shape[0]
	for i in xrange(n):
		i_topic = cl_idx[i]
		for j in xrange(n):
			if i==j:
				continue
			elif i_topic==cl_idx[j]:
				distanceMatrix[i,j]*=sameTopicWeight
			else:
				distanceMatrix[i,j]*=differentTopicWeight
		distanceMatrix[i,:] /= distanceMatrix[i,:].max()

	return np.round(distanceMatrix,decimals=4)


# Routing
# @app.before_first_request
def before__first_request_():
	global eng
	global mappedX
	global cl_idx
	global Wtopk
	global voca
	global distanceMatrix
	global distanceMatrix_main


	tic = time.time()
	print 'Starting matlab - ',
	eng = engine.start_matlab()
	eng.cd(os.path.dirname(os.getcwd()))
	print "%.4f" % (time.time()-tic)

	tic = time.time()
	print "Get data - ",
	[mappedX, cl_idx, Wtopk_idx,voca] = eng.main_topic(nargout=4)
	distanceMatrix = io.loadmat('./../tdm2.mat')['DD']
	print "%.4f" % (time.time()-tic)


	tic = time.time()
	print "Calculate data - ",

	Wtopk = []
	for idxArray in Wtopk_idx:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk.append(tempArray)

	cl_idx = cl_idx[0]


	cl_idx = cl_idx[0:1000]
	distanceMatrix = distanceMatrix[0:1000,0:1000]


	sameTopicWeight = 0.8
	differentTopicWeight = 1.15
	distanceMatrix_main = supervisedTSNE(distanceMatrix, cl_idx,
		sameTopicWeight=sameTopicWeight, differentTopicWeight=differentTopicWeight)

	print "%.4f" % (time.time()-tic)

	distanceMatrix_main = distanceMatrix_main.tolist()

@app.teardown_request
def teardown_request(exception):
	print('Teardown arose!'.format(exception))


@app.route('/get_subTopic')
def get_subTopic():
	global eng
	global voca
	global distanceMatrix

	idx = json.loads(request.args.get('idx'))

	[mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = eng.sub_topic(idx,nargout=3)

	Wtopk_sub = []
	for idxArray in Wtopk_idx_sub:
		tempArray = []
		for topicIdx in idxArray:
			tempArray.append(voca[int(topicIdx)-1])
		Wtopk_sub.append(tempArray)

	cl_idx_sub = cl_idx_sub[0]
	cl_idx_sub = np.array(cl_idx_sub).tolist()


	distanceMatrix_sub = distanceMatrix[idx,:][:,idx]

	sameTopicWeight = 0.8
	differentTopicWeight = 1.15
	distanceMatrix_sub = supervisedTSNE(distanceMatrix_sub, cl_idx,
		sameTopicWeight=sameTopicWeight, differentTopicWeight=differentTopicWeight)
	distanceMatrix_sub = distanceMatrix_sub.tolist()



	return json.dumps({'distanceMatrix':distanceMatrix_sub, 'cl_idx_sub':cl_idx_sub, 'Wtopk_sub':Wtopk_sub})


@app.route('/get_subTopic_tsne')
def get_subTopic_tsne():
	global eng
	global voca
	idx = json.loads(request.args.get('idx'))

	[mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = eng.sub_topic_tsne(idx,nargout=3)

	print mappedX_sub

	Wtopk_sub = []
	for idxArray in Wtopk_idx_sub:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk_sub.append(tempArray)

	cl_idx_sub = cl_idx_sub[0]

	mappedX_sub = np.array(mappedX_sub).tolist()
	cl_idx_sub = np.array(cl_idx_sub).tolist()

	return json.dumps({'mappedX_sub':mappedX_sub, 'cl_idx_sub':cl_idx_sub, 'Wtopk_sub':Wtopk_sub})

# keyword 입력받음
@app.route('/')
def form():
	global cl_idx
	global Wtopk
	global distanceMatrix_main

	return render_template('tsne.html', cl_idx=cl_idx, Wtopk= Wtopk, distanceMatrix=distanceMatrix_main)
 

# Execute the main program
if __name__ == '__main__':
	before__first_request_()
	app.run(host='0.0.0.0',port=5004)
