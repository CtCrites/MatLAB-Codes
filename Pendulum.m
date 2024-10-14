   %   P   e   n   d   u   l       -       P   r   o   g   r   a   m       t   o       c   o   m   p   u   t   e       t   h   e       m   o   t   i   o   n       o   f       a       s   i   m   p   l   e       p   e   n   d   u   l   u   m      
   %   u   s   i   n   g       E   u   l   e   r   '   s       M   e   t   h   o   d      
   c   l   e   a   r       a   l   l   ;       %   C   l   e   a   r       m   e   m   o   r   y       a   n   d       p   r   i   n   t       t   h   e       h   e   a   d   e   r      
      
      
   %   S   e   t       i   n   i   t   i   a   l       v   a   l   u   e   s      
   %   S   e   t       c   o   n   s   t   a   n   t   s       a   n   d       d   e   f   i   n   e       o   t   h   e   r       v   a   r   i   a   b   l   e   s      
   g   _   o   v   e   r   _   L       =       1   ;      
      
   t   a   u       =       0   .   0   0   0   8   ;                                                                                                   %   S   i   z   e       o   f       t   i   m   e       s   t   e   p      
   n   s   t   e   p       =       3   0   0   0   0   ;                                                                                               %   N   u   m   b   e   r       o   f       s   t   e   p   s      
   t   i   m   e       =       0   ;                                                                                                               %   S   e   t   s       s   t   a   r   t   i   n   g       t   i   m   e      
   i   r   e   v       =       0   ;                                                                                                               %   S   e   t   s       s   t   a   r   t   i   n   g       r   e   v   o   l   u   t   i   o   n      
      
   %   P   l   o   t       1      
   t   h   e   t   a   0       =       1   5   ;                                                                                                   %   i   n   i   t   i   a   l       a   n   g   l   e      
   t   h   e   t   a       =       t   h   e   t   a   0   *   p   i   /   1   8   0   ;                                                           %   C   o   n   v   e   r   t       i   n   i   t   i   a   l       a   n   g   l   e       t   o       r   a   d   i   a   n   s      
   o   m   e   g   a       =       0   .   0   1   ;                                                                                               %   I   n   i   t   i   a   l       t   a   n   g   e   n   t   i   a   l       v   e   l   o   c   i   t   y      
   a   c   c   e   l       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   )   ;       %   S   e   t   s       a   c   c   e   l   e   r   a   t   i   o   n      
      
   %   P   l   o   t       2      
   t   h   e   t   a   0   2       =       6   5   ;      
   t   h   e   t   a   2       =       t   h   e   t   a   0   2   *   p   i   /   1   8   0   ;      
   o   m   e   g   a   2       =       0   .   0   1   ;      
   a   c   c   e   l   2       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   2   )   ;      
      
   %   P   l   o   t       3      
   t   h   e   t   a   0   3       =       9   0   ;      
   t   h   e   t   a   3       =       t   h   e   t   a   0   3   *   p   i   /   1   8   0   ;      
   o   m   e   g   a   3       =       1   .   4   2   3   ;      
   a   c   c   e   l   3       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   3   )   ;      
      
   %   P   l   o   t       4      
   t   h   e   t   a   0   4       =       9   0   ;      
   t   h   e   t   a   4       =       t   h   e   t   a   0   4   *   p   i   /   1   8   0   ;      
   o   m   e   g   a   4       =       1   .   4   2   3   ;      
   a   c   c   e   l   4       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   4   )   ;      
      
   %   S   e   t       u   p       l   o   o   p      
   f   o   r       i   s   t   e   p       =       1   :   n   s   t   e   p      
                   t   _   p   l   o   t   (   i   s   t   e   p   )       =       t   i   m   e   ;      
                   t   h   _   p   l   o   t   (   i   s   t   e   p   )       =       t   h   e   t   a   *   1   8   0   /   p   i   ;      
                   o   m   _   p   l   o   t   (   i   s   t   e   p   )       =       o   m   e   g   a       +       t   a   u   *   a   c   c   e   l   ;      
                   t   i   m   e       =       t   i   m   e       +       t   a   u   ;      
      
                   t   h   _   p   l   o   t   2   (   i   s   t   e   p   )       =       t   h   e   t   a   2   *   1   8   0   /   p   i   ;      
                   o   m   _   p   l   o   t   2   (   i   s   t   e   p   )       =       o   m   e   g   a   2       +       t   a   u   *   a   c   c   e   l   2   ;      
      
                   t   h   _   p   l   o   t   3   (   i   s   t   e   p   )       =       t   h   e   t   a   3   *   1   8   0   /   p   i   ;      
                   o   m   _   p   l   o   t   3   (   i   s   t   e   p   )       =       o   m   e   g   a   3       +       t   a   u   *   a   c   c   e   l   3   ;      
      
                   t   h   _   p   l   o   t   4   (   i   s   t   e   p   )       =       t   h   e   t   a   4   *   1   8   0   /   p   i   ;      
                   o   m   _   p   l   o   t   4   (   i   s   t   e   p   )       =       o   m   e   g   a   4       +       t   a   u   *   a   c   c   e   l   4   ;      
      
                   t   h   e   t   a       =       t   h   e   t   a       +       t   a   u   *   o   m   e   g   a   ;      
                   o   m   e   g   a       =       o   m   e   g   a       +       t   a   u   *   a   c   c   e   l   ;      
                   a   c   c   e   l       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   )   ;      
      
                   t   h   e   t   a   2       =       t   h   e   t   a   2       +       t   a   u   *   o   m   e   g   a   2   ;      
                   o   m   e   g   a   2       =       o   m   e   g   a   2       +       t   a   u   *   a   c   c   e   l   2   ;      
                   a   c   c   e   l   2       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   2   )   ;      
      
                   t   h   e   t   a   3       =       t   h   e   t   a   3       +       t   a   u   *   o   m   e   g   a   3   ;      
                   o   m   e   g   a   3       =       o   m   e   g   a   3       +       t   a   u   *   a   c   c   e   l   3   ;      
                   a   c   c   e   l   3       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   3   )   ;      
      
                   t   h   e   t   a   4       =       t   h   e   t   a   4       +       t   a   u   *   o   m   e   g   a   4   ;      
                   o   m   e   g   a   4       =       o   m   e   g   a   4       +       t   a   u   *   a   c   c   e   l   4   ;      
                   a   c   c   e   l   4       =       -   g   _   o   v   e   r   _   L   *   s   i   n   (   t   h   e   t   a   4   )   ;                      
                      
      
   e   n   d      
      
   f   i   g   u   r   e   (   g   c   f   )   ;      
   p   l   o   t   (   t   h   _   p   l   o   t   ,       o   m   _   p   l   o   t   ,       '   .   '   ,       t   h   _   p   l   o   t   2   ,       o   m   _   p   l   o   t   2   ,       '   -   '   ,       t   h   _   p   l   o   t   3   ,       o   m   _   p   l   o   t   3   ,       '   -   -   '   ,       t   h   _   p   l   o   t   4   ,       o   m   _   p   l   o   t   4   ,       '   -   '   )      
   y   l   a   b   e   l   (   "   \   o   m   e   g   a       (   m   /   s   )   "   )   ;       x   l   a   b   e   l   (   "   \   t   h   e   t   a       (   d   e   g   r   e   e   s   )   "   )   ;      
   t   i   t   l   e   (   "   \   t   h   e   t   a       v   .       \   o   m   e   g   a   "   )   ;      
      
      
      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ��