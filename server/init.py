# -*- coding: utf-8 -*-

from werkzeug.contrib.profiler import ProfilerMiddleware
from flask import Flask, request, g, render_template
from flask.ext.triangle import Triangle
from flask.ext.socketio import SocketIO, emit
from scipy import sparse, io
import numpy as np
from matlab import engine
import os, json, time
from flask.ext.cors import CORS


# Configuration
app = Flask(__name__, static_path='/static')	
Triangle(app)
CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

app.config['DEBUG'] = True
app.config.from_object(__name__)


app.config['PROFILE'] = True
app.wsgi_app = ProfilerMiddleware(app.wsgi_app, restrictions=[10])
socketio = SocketIO(app)
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

	return np.round(distanceMatrix,decimals=3)


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
	eng.main_topic2(nargout=0)

	mappedX = eng.workspace['mappedX']
	cl_idx = eng.workspace['cl_idx']
	Wtopk_idx = eng.workspace['Wtopk_idx']
	voca = eng.workspace['dict']
	print "%.4f" % (time.time()-tic)
	
	distanceMatrix = io.loadmat('./../newtdm2.mat')['DD']
	
	tic = time.time()
	print "Calculate data - ",

	Wtopk = []
	for idxArray in Wtopk_idx:
		tempArray = []
		for idx in idxArray:
			tempArray.append(voca[int(idx)-1])
		Wtopk.append(tempArray)

	cl_idx = cl_idx[0]


	# cl_idx = cl_idx
	# distanceMatrix = distanceMatrix

	# cl_idx = cl_idx
	# distanceMatrix = distanceMatrix

	sameTopicWeight = 0.9
	differentTopicWeight = 1.1
	distanceMatrix_main = supervisedTSNE(distanceMatrix, cl_idx,
		sameTopicWeight=sameTopicWeight, differentTopicWeight=differentTopicWeight)

	print "%.4f" % (time.time()-tic)

	distanceMatrix_main = distanceMatrix_main.tolist()

@app.teardown_request
def teardown_request(exception):
	pass


@app.route('/get_subTopic')
def get_subTopic():
	global eng
	global voca
	global distanceMatrix

	idx = json.loads(request.args.get('idx'))

	eng.workspace['idx'] = idx
	eng.sub_topic(nargout=0)
	mappedXP_sub = eng.workspace['mappedX_sub']
	cl_idx_sub = eng.workspace['cl_idx_sub']
	Wtopk_idx_sub = eng.workspace['Wtopk_idx_sub']
	k_sub = eng.workspace['k_sub'] # number of topics that will be shown
	

	idx = [i-1 for i in idx]

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
	differentTopicWeight = 1.2

	distanceMatrix_sub = supervisedTSNE(distanceMatrix_sub, cl_idx_sub,
		sameTopicWeight=sameTopicWeight, differentTopicWeight=differentTopicWeight)
	distanceMatrix_sub = distanceMatrix_sub.tolist()



	return json.dumps({'distanceMatrix':distanceMatrix_sub, 'cl_idx_sub':cl_idx_sub, 'Wtopk_sub':Wtopk_sub})



###### socket test code
@socketio.on('connect', namespace='/subtopic')
def connect():
	print "connected"

@socketio.on('disconnect', namespace='/subtopic')
def disconnect():
	print "disconnected"

@socketio.on('get_subTopic', namespace='/subtopic')
def get_subTopic_(message):
	global eng
	global voca
	global distanceMatrix
	
	idx = message['idx']
	socketId = message['socketId']
	print "get Request - %s" % socketId
	sameTopicWeight = 0.75
	differentTopicWeight = 1.25
	
	eng.workspace['idx'] = idx
	print "before subtopic"
	eng.sub_topic(nargout=0)
	print "after subtopic"
	k_sub = int(eng.workspace['k_sub'])
	sub_k = int(eng.workspace['sub_k'])	
	idx = [i-1 for i in idx]
	distanceMatrix_sub = distanceMatrix[idx,:][:,idx]
	print "before iteration"
	
	if k_sub-sub_k==0:
		iterNum = 1
	else:
		iterNum = k_sub-sub_k

	for i in xrange(1,iterNum+1):
		print i
		eng.workspace['i'] = i
		eng.sub_topic_ith_Iter(nargout=0)

		cl_idx_sub = eng.workspace['cl_idx_sub']
		Wtopk_idx_sub = eng.workspace['Wtopk_idx_sub']

		Wtopk_sub = []
		for idxArray in Wtopk_idx_sub:
			tempArray = []
			for topicIdx in idxArray:
				tempArray.append(voca[int(topicIdx)-1])
			Wtopk_sub.append(tempArray)

		cl_idx_sub = cl_idx_sub[0]
		cl_idx_sub = np.array(cl_idx_sub).tolist()

		distanceMatrix_sub_ = supervisedTSNE(distanceMatrix_sub, cl_idx_sub,
			sameTopicWeight=sameTopicWeight, differentTopicWeight=differentTopicWeight)
		distanceMatrix_sub_ = distanceMatrix_sub.tolist()

		emit('result data'+socketId, {'distanceMatrix':distanceMatrix_sub_, 'cl_idx_sub':cl_idx_sub, 'Wtopk_sub':Wtopk_sub})
	# time.sleep(15)
	return 1




#########################
# not used in gather plot
@app.route('/get_subTopic_tsne')
def get_subTopic_tsne():
	global eng
	global voca
	idx = json.loads(request.args.get('idx'))

	eng.sub_topic(nargout=0)
	mappedXP_sub = eng.workspace['mappedX_sub']
	cl_idx_sub = eng.workspace['cl_idx_sub']
	Wtopk_idx_sub = eng.workspace['Wtopk_idx_sub']
	k_sub = eng.workspace['k_sub'] # number of topics that will be shown

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
	socketio.run(app,host='0.0.0.0',port=5004, debug=True)
