  I)  ^   k820309    )          19.0        .<_                                                                                                          
       pes.f90 NNPES                      @                              
                         @               @                '             $      #RX    #M    #L    #MAXS    #RY    #S    #W 	   #B 
   #WC    #BC    #Y    #DYDX    #A    #N    #DYDN    #DYDP    #DNDX    #DADX    #DXDYDN    #DPDYDX    #ACTIVATION    #STRUCTFL    #INIT    #OUTPUT    #WB_INIT     #WB_COPY %   #WB_REC (   #WB +   #WB_UPDATE 0   #CAL_DYDN 4   #CAL_DYDP 7   #CAL_DADX :   #CAL_DXDYDN =   #CAL_DPDYDX @   #SAVENN C   #CLEAN H                �                                                               �                                                              �                                                              �                                                              �                                                             �                                           (                             &                                                      �                              	            p                 
            &                   &                   &                                                      �                              
            �                 
            &                   &                                                      �                                          H             	   
            &                   &                   &                                                      �                                          �             
   
            &                   &                                                      �                                                           
            &                                                      �                                          h                
            &                   &                                                      �                                          �                
            &                   &                                                      �                                          (                
            &                   &                                                      �                                          �                
            &                   &                   &                                                      �                                                           
            &                   &                                                      �                                          `                
            &                   &                   &                                                      �                                          �                
            &                   &                   &                                                      �                                          P                
            &                   &                   &                   &                                                      �                                          �                
            &                   &                   &                                           .           �                                          X      c                      &                                                                �                                   c       �             1         �   � $                      �                        #NN_INIT    #         @     @                                                #THIS                                                                 #NN    1         �   � $                      �                        #NN_OUTPUT    #         @     @                                                #THIS    #X                                                                 #NN             
                                                     
    p          5 8 O#NN     p        U            5 8 O#NN     p        U                          1         �   � $                      �                         #INIT_WB !   #         @     @                            !                    #THIS "   #X1 #   #X2 $                                             "                   #NN              
                                 #     
                
                                 $     
      1         �   � $                      �      %                  #COPY_WB &   #         @     @                            &                    #THIS '                                             '                   #NN    1         �   � $                      �      (                  #REC_WB )   #         @     @                            )                    #THIS *                                             *                   #NN    1         �   � $                      �      +                  #WB_INOUT ,   #         @     @                            ,                    #THIS -   #ID .   #DX /                                             -                   #NN              
                                  .                                                    /                    
     p          5 8 O#NN     p        U            5 8 O#NN     p        U                          1         �   � $                      �      0                  #UPDATE_WB 1   #         @     @                            1                    #THIS 2   #DX 3                                             2                   #NN             
                                 3                    
    p          5 8 O#NN     p        U            5 8 O#NN     p        U                          1         �   � $                      �      4                  #CALC_DYDN 5   #         @     @                            5                    #THIS 6                                             6                   #NN    1         �   � $                      �      7              	    #CALC_DYDP 8   #         @     @                            8                    #THIS 9                                             9                   #NN    1         �   � $                      �      :               
    #CALC_DADX ;   #         @     @                            ;                    #THIS <                                             <                   #NN    1         �   � $                      �      =             !     #CALC_DXDYDN >   #         @     @                            >                    #THIS ?                                             ?                   #NN    1         �   � $                      �      @             "     #CALC_DPDYDX A   #         @     @                            A                    #THIS B                                             B                   #NN    1         �   � $                      �      C             #     #SAVEPARA D   #         @     @                            D                    #THIS E   #FILENAME F   #IO G                                             E                   #NN              
                                F                    1           
                                  G           1         �   � $                      �      H             $     #CLEAN_NN I   #         @     @                            I                    #THIS J                                             J                   #NN                                                K                                                        L           #NN                                                M                       @                                 N                              @                           O     '                   #CTYPE P   #NPERM Q   #ALIST R   #ACOEF S   #SCALING T   #SCCOEF U                �                               P                                �                               Q                               �                               R                               p          p          p            p          p                                       �                               S            �                   p          p            p                                       �                               T                               �                               U                            
  p          p            p                                            @                           V     'H                   #NC W   #TERMS X                �                               W                                �                               X     (                          p          p          p 
           p          p 
                                  @ @                               Y                                   &                                           #CDEF O            @ @                               Z            H                       &                   &                                           #HDFORM V              @                                 [                       @                                 \            #         @                                   ]                        �         fn#fn    �   @   J   NN_CLASS    �   �      NN+NN_CLASS    �  H   a   NN%RX+NN_CLASS      H   a   NN%M+NN_CLASS    T  H   a   NN%L+NN_CLASS !   �  H   a   NN%MAXS+NN_CLASS    �  H   a   NN%RY+NN_CLASS    ,  �   a   NN%S+NN_CLASS    �  �   a   NN%W+NN_CLASS    �  �   a   NN%B+NN_CLASS    0  �   a   NN%WC+NN_CLASS    �  �   a   NN%BC+NN_CLASS    �  �   a   NN%Y+NN_CLASS !   4  �   a   NN%DYDX+NN_CLASS    �  �   a   NN%A+NN_CLASS    �	  �   a   NN%N+NN_CLASS !   8
  �   a   NN%DYDN+NN_CLASS !   �
  �   a   NN%DYDP+NN_CLASS !   �  �   a   NN%DNDX+NN_CLASS !   l  �   a   NN%DADX+NN_CLASS #   0  �   a   NN%DXDYDN+NN_CLASS #     �   a   NN%DPDYDX+NN_CLASS '   �  �   a   NN%ACTIVATION+NN_CLASS %   l  P   a   NN%STRUCTFL+NN_CLASS !   �  U   a   NN%INIT+NN_CLASS !     R      NN_INIT+NN_CLASS &   c  P   a   NN_INIT%THIS+NN_CLASS #   �  W   a   NN%OUTPUT+NN_CLASS #   
  Y      NN_OUTPUT+NN_CLASS (   c  P   a   NN_OUTPUT%THIS+NN_CLASS %   �  �   a   NN_OUTPUT%X+NN_CLASS $   �  U   a   NN%WB_INIT+NN_CLASS !   �  b      INIT_WB+NN_CLASS &   B  P   a   INIT_WB%THIS+NN_CLASS $   �  @   a   INIT_WB%X1+NN_CLASS $   �  @   a   INIT_WB%X2+NN_CLASS $     U   a   NN%WB_COPY+NN_CLASS !   g  R      COPY_WB+NN_CLASS &   �  P   a   COPY_WB%THIS+NN_CLASS #   	  T   a   NN%WB_REC+NN_CLASS     ]  R      REC_WB+NN_CLASS %   �  P   a   REC_WB%THIS+NN_CLASS    �  V   a   NN%WB+NN_CLASS "   U  b      WB_INOUT+NN_CLASS '   �  P   a   WB_INOUT%THIS+NN_CLASS %     @   a   WB_INOUT%ID+NN_CLASS %   G  �   a   WB_INOUT%DX+NN_CLASS &     W   a   NN%WB_UPDATE+NN_CLASS #   v  Z      UPDATE_WB+NN_CLASS (   �  P   a   UPDATE_WB%THIS+NN_CLASS &      �   a   UPDATE_WB%DX+NN_CLASS %   �  W   a   NN%CAL_DYDN+NN_CLASS #   O  R      CALC_DYDN+NN_CLASS (   �  P   a   CALC_DYDN%THIS+NN_CLASS %   �  W   a   NN%CAL_DYDP+NN_CLASS #   H  R      CALC_DYDP+NN_CLASS (   �  P   a   CALC_DYDP%THIS+NN_CLASS %   �  W   a   NN%CAL_DADX+NN_CLASS #   A  R      CALC_DADX+NN_CLASS (   �  P   a   CALC_DADX%THIS+NN_CLASS '   �  Y   a   NN%CAL_DXDYDN+NN_CLASS %   <  R      CALC_DXDYDN+NN_CLASS *   �  P   a   CALC_DXDYDN%THIS+NN_CLASS '   �  Y   a   NN%CAL_DPDYDX+NN_CLASS %   7  R      CALC_DPDYDX+NN_CLASS *   �  P   a   CALC_DPDYDX%THIS+NN_CLASS #   �  V   a   NN%SAVENN+NN_CLASS "   /  h      SAVEPARA+NN_CLASS '   �  P   a   SAVEPARA%THIS+NN_CLASS +   �  L   a   SAVEPARA%FILENAME+NN_CLASS %   3   @   a   SAVEPARA%IO+NN_CLASS "   s   V   a   NN%CLEAN+NN_CLASS "   �   R      CLEAN_NN+NN_CLASS '   !  P   a   CLEAN_NN%THIS+NN_CLASS    k!  @       NATOMS    �!  H       ANN    �!  @       NCOORD    3"  @       NSTATES    s"  �       CDEF    #  H   a   CDEF%CTYPE    P#  H   a   CDEF%NPERM    �#  �   a   CDEF%ALIST    T$  �   a   CDEF%ACOEF    �$  H   a   CDEF%SCALING    8%  �   a   CDEF%SCCOEF    �%  c       HDFORM    7&  H   a   HDFORM%NC    &  �   a   HDFORM%TERMS    ;'  �       COORD    �'  �       HDDEF    �(  @       LASTLAYER    �(  @       NOOPC    )  H       READCOORDINPUT 