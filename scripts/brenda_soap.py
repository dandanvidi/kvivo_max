# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 17:14:46 2015

@author: noore
"""
from SOAPpy import SOAPProxy
import hashlib

class BrendaSoap(object):
    
    def __init__(self, username, password):
        self.client = self._GetClient()
        self.My_Username = username
        password = hashlib.sha256(password).hexdigest()
        self.account = username+","+password+","
                    
    def _GetClient(use_soap=True):
        if use_soap:
            # create the SOAP client object
            #endpointURL = "http://www.brenda-enzymes.org/soap2/brenda_server.php"
            endpointURL = "http://www.brenda-enzymes.org/soap/brenda_server.php"
            return SOAPProxy(endpointURL)
        else:
            import WSDL
            # create the WSDL SOAP client object
            wsdl = "http://www.brenda-enzymes.org/soap2/brenda.wsdl"
            return WSDL.Proxy(wsdl)

    @staticmethod
    def _ParseTableResult(res):
        row_dicts = []
        print res
        if not res:
            return row_dicts
        for entry in res.split('!'):
            if not entry:
                continue
            if entry[-1] == '#':
                entry = entry[:-1]
            row_dict = dict(map(lambda s:s.split('*', 1), entry.split('#')))
            row_dicts.append(row_dict)
        return row_dicts
    
    @staticmethod
    def _ParseListResult(res):
        return res.split('!')
        
    @staticmethod
    def _MakeQuery(ec=None, organism=None):
        s = []
        if ec is not None:
            s += ['ecNumber*' + ec]
        if organism is not None:
            s += ['organism*' + organism]
        if s == []:
            raise ValueError('either EC number of Organism must be provided')
        return '#'.join(s)
        
    def GetTurnoverNumber(self, ec=None, organism=None):
        parameters = self.account+BrendaSoap._MakeQuery(ec, organism)
        res = self.client.getTurnoverNumber(parameters)
        return BrendaSoap._ParseTableResult(res)
        
    def GetActivatingCompound(self, ec=None, organism=None):
        parameters = self.account+BrendaSoap._MakeQuery(ec, organism)
        res = self.client.getActivatingCompound(parameters)
        return BrendaSoap._ParseTableResult(res)
        
    def GetInhibitors(self, ec=None, organism=None):
        parameters = self.account+BrendaSoap._MakeQuery(ec, organism)
        res = self.client.getInhibitors(parameters)
        return BrendaSoap._ParseTableResult(res)

    def GetEcNumbersFromActivatingCompound(self):
        res = self.client.getEcNumbersFromActivatingCompound()
        return BrendaSoap._ParseListResult(res)
        
    def GetOrganismsFromActivatingCompound(self):
        res = self.client.getOrganismsFromActivatingCompound()
        return BrendaSoap._ParseListResult(res)
        
    def GetEcNumbersFromInhibitors(self):
        res = self.client.getEcNumbersFromInhibitors()
        return BrendaSoap._ParseListResult(res)
        
    def GetOrganismsFromInhibitors(self):
        res = self.client.getOrganismsFromInhibitors()
        return BrendaSoap._ParseListResult(res)
    
