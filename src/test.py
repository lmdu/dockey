from collections import OrderedDict

class AttrDict(OrderedDict):
    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, attr, val):
        self[attr] = val

o = AttrDict({
    'c': 1,
    'd': 2,
    'e': 3
})

for k in o:
    print(k)
