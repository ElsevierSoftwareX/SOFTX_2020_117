from CodernityDB.tree_index import TreeBasedIndex

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