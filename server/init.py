# -*- coding: utf-8 -*-

from werkzeug.contrib.profiler import ProfilerMiddleware
from flask import Flask, request, g, render_template
from flask.ext.triangle import Triangle
from scipy import sparse, io
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
from matlab import engine
import os, json
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


# Routing
# @app.before_first_request
def before__first_request():
	global eng
	global mappedX
	global cl_idx
	global Wtopk
	global voca
	global distanceMatrix


	print 'Starting matlab'
	eng = engine.start_matlab()
	eng.cd(os.path.dirname(os.getcwd()))
	[mappedX, cl_idx, Wtopk_idx,voca] = eng.main_topic(nargout=4)

	distanceMatrix = io.loadmat('./../tdm2.mat')['DD']

	Wtopk = []
	for idxArray in Wtopk_idx:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk.append(tempArray)

	cl_idx = cl_idx[0]


	sameTopicWeight = 0.96
	differentTopicWeight = 1.04

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

	distanceMatrix = np.round(distanceMatrix,decimals=4)

	print 'Server Ready'
	

	cl_idx = cl_idx[0:1000]
	distanceMatrix = distanceMatrix[0:1000,0:1000]

	distanceMatrix = distanceMatrix.tolist()

@app.teardown_request
def teardown_request(exception):
	print('Teardown arose!'.format(exception))


@app.route('/get_subTopic')
def get_subTopic():
	global eng
	global voca
	idx = json.loads(request.args.get('idx'))

	[mappedX_sub, cl_idx_sub, Wtopk_idx_sub] = eng.sub_topic(idx,nargout=3)

	print mappedX_sub

	Wtopk_sub = []
	for idxArray in Wtopk_idx_sub:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk_sub.append(tempArray)

	cl_idx_sub = cl_idx_sub[0]

	# mappedX_sub = np.array(mappedX_sub).tolist()
	cl_idx_sub = np.array(cl_idx_sub).tolist()

	return json.dumps({'mappedX_sub':mappedX_sub, 'cl_idx_sub':cl_idx_sub, 'Wtopk_sub':Wtopk_sub})


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
	global mappedX
	global cl_idx
	global Wtopk
	global distanceMatrix

	return render_template('tsne.html', cl_idx=cl_idx, Wtopk= Wtopk, distanceMatrix=distanceMatrix)
 
before__first_request()

# Execute the main program
if __name__ == '__main__':
	app.run(host='0.0.0.0',port=5004)
