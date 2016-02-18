# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 13:14:45 2016

@author: minjeong1642
"""
def t_sne():
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.cd('D:')
    eng.cd('D:\\Minjeong\\nmf_topicmodeling')
    [mappedX, mappedX_sub, cl_idx, cl_idx_sub, Wtopk_idx, Wtopk_idx_Sub]=eng.script_runme_rank2_python(nargout=6)
    return mappedX, mappedX_sub, cl_idx, cl_idx_sub, Wtopk_idx, Wtopk_idx_Sub
