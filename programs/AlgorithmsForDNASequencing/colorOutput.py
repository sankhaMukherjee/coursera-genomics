class Color:

    def __init__(self, s):
        self.s = s

    def warn(self):
        return '\x1b[0;31;40m'+self.s+'\x1b[0m'

def cStr(s1, s2):
    '''colorize s1 based on s2
    
    giveb 2 string, this will give string s1 colorized
    with red where this is not the same as s2. If s1
    and s2 are not the same length, a '\xb0' wlill be added
    to the end
    
    Arguments:
        s1 {str} -- string for comparison
        s2 {str} -- string for comparison
    '''

    remain = ''
    if len(s1) < len(s2):
        remain = Color('-'*( len(s2) - len(s1) )).warn()
        
    if len(s2) < len(s1):
        remain = Color(s1[len(s2):]).warn()

    result = [ (a if a == b else Color(a).warn()) for a, b in zip(s1, s2)]
    result = ''.join(result)
    result += remain

    return result

if __name__ == '__main__':
    
    print(cStr('sankha', 'saxkxa'))
    print(cStr('sankhaSubhraMukherjee', 'sankha'))
    print(cStr('sankha', 'sankhaSubhraMukherjee'))
