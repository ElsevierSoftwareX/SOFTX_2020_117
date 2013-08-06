from CodernityDB.tree_index import TreeBasedIndex
from CodernityDB.hash_index import HashIndex
from hashlib import md5

class TestIDIndex(TreeBasedIndex):

    def __init__(self, *args, **kwargs):
        kwargs['key_format'] = 'I'
        super(TestIDIndex, self).__init__(*args, **kwargs)

    def make_key_value(self, data):
        if data['t'] == 'test':
            key = int(data['test_id'])
            return key, None
        else:
            return None
        
    def make_key(self, key):
        return key

class SrcCommitIndex(HashIndex):

    def __init__(self, *args, **kwargs):
        kwargs['key_format'] = '40s'
        super(SrcCommitIndex, self).__init__(*args, **kwargs)

    def make_key_value(self, data):
        if data['t'] == 'test':
            key = str(data['git_commit'])
            return key, None
        else:
            return None
        
    def make_key(self, key):
        return key
        
        
class KatTestIndex(HashIndex):

    def __init__(self, *args, **kwargs):
        kwargs['key_format'] = '32s'
        super(KatTestIndex, self).__init__(*args, **kwargs)

    def make_key_value(self, data):
        if data['t'] == 'kattest':
            key = str(data['suite'])  + "_" + str(data['kat'])
            return md5(key).hexdigest(), None
        else:
            return None

    def make_key(self, key):
        return md5(key).hexdigest()
        