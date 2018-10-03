#!/usr/bin/env python

# Canvas module for simple drawing
# Quintin Cutts
# 23 - 8 - 06

from Tkinter import *
import threading
import time
import exceptions

class WindowGone(exceptions.Exception):
    def __init__(self, args=[]):
        self.args = args

class RawCanvas:
    mainLoopRunning = False
    def __init__(self, canvas):
        self._canvas = canvas
    def create_rectangle( self, x1, y1, x2, y2, *kw ):
        r = self._canvas.create_rectangle( x1, y1, x2, y2, kw )
        self._canvas._root().update()
        return r
    def create_arc( self, x1, y1, x2, y2, *kw ):
        r = self._canvas.create_arc( x1, y1, x2, y2, kw )
        self._canvas._root().update()
        return r        
    def create_line( self, x1, y1, x2, y2, *kw ):
        r = self._canvas.create_line( x1, y1, x2, y2, kw)
        self._canvas._root().update() 
        return r       
    def create_oval( self, x1, y1, x2, y2, *kw ):
        r = self._canvas.create_oval( x1, y1, x2, y2, kw )
        self._canvas._root().update()
        return r        
    def create_text( self, x1, y1, *kw ):
        r = self._canvas.create_text( x1, y1, kw )
        self._canvas._root().update() 
        return r       
    def move( self, tagOrId, xInc, yInc ):
        self._canvas.move( tagOrId, xInc, yInc )
        self._canvas._root().update()        
    def delete( self, tagOrId ):
        self._canvas.delete( tagOrId )
        self._canvas._root().update()        
         
    def complete( self ):
        global _can
        print _can
        _can.destroy()
        print "after can.destroy()"
##        print "in complete"
##        print self.mainLoopRunning
##        self._canvas.unbind("<Button-1>")
##        self._canvas.bind("<Button-1>", _can.destroy)
##        wait( 0.5 )
##        print "in between"
##        self._canvas.create_line( 10,10,100,100 )
##        print "after draw line"
##        x = self._canvas._root()
##        print " after root call"
##        y = x.title
##        print "after getting title"
##        print y()
##        print "after printing title"
##        z = y( "new title" )
##        x.update()
##        #self._canvas._root().title( "Click mouse to end" )
##        print "before self.run"
##        self.run()
            
    def run( self ):
        if not self.mainLoopRunning:
            self.mainLoopRunning = True
            try:
                self._canvas._root().mainloop()
            except WindowGone:
                pass
         

_root = None
_canvas = None      # This is the real Python Tkinter canvas
_can = None         # This is the Glasgow canvas
_myThreads = []
_hadCan = False
_blockCalls = False

class Can( RawCanvas ):
    def __init__( self ):
        global _root, _canvas
        if _root is None:
            _root = Tk()
        if _canvas is None:
            _canvas = Canvas( _root, background = "white" )
            _canvas.pack(expand=1, fill="both" )
        RawCanvas.__init__( self, _canvas )

        _root.iconify()
        _root.update()
        _root.deiconify()
        #_root.lift()

        _root.update()
        
        def onWinClose():
            global _blockCalls
            _blockCalls = True
            time.sleep(0.5 )
            _root.destroy()
        
        
        _root.protocol("WM_DELETE_WINDOW",self.destroy )

    def destroy( event, extra=None ):
        global _blockCalls, _root
        _blockCalls = True
        time.sleep( 0.5 )
        _root.destroy()

def _getCanvas():
    global _can, _hadCan, _blockCalls
    can = _can
    if (_hadCan and not can) or _blockCalls:
        raise WindowGone
    if not can:
        _can = can = Can()
        _hadCan = True
    return can

no_current_keyhandler_call = True    # Concurrency control - stops multiple simultaneous calls of the handler
no_current_mousehandler_call = True
##########################################################
# These are the only visible functions out of the module

def create_rectangle( x1, y1, x2, y2, **kw ):
    return _getCanvas().create_rectangle( x1, y1, x2, y2, kw )
def create_arc( x1, y1, x2, y2, **kw ):
    return _getCanvas().create_arc( x1, y1, x2, y2, kw )
def create_line( x1, y1, x2, y2, **kw ):
    return _getCanvas().create_line( x1, y1, x2, y2, kw )
def create_oval( x1, y1, x2, y2, **kw ):
    return _getCanvas().create_oval( x1, y1, x2, y2, kw )
def create_text( x1, y1, **kw ):
    return _getCanvas().create_text( x1, y1, kw )
def move( tagOrId, xInc, yInc ):
    _getCanvas().move( tagOrId, xInc, yInc )
def wait( t1 ):
    time.sleep( t1 )
def delete( tagOrId ):
    _getCanvas().delete( tagOrId )
def set_size( x, y ):
    _getCanvas()._canvas.config( width = x, height = y )
def complete():
    _getCanvas().complete()
def run():
    _getCanvas().run()
#def quitCanvas():
def runGraphicsFn( g ):
    global _myThreads
    if _myThreads == []:
        create_rectangle( 1,1,2,2,outline="black" )  #ensure a canvas has been created

    def gWrap():
        try:
            print "starting thread"
            g()
        except WindowGone:
            pass

    newThread = threading.Thread( target = gWrap )
    _myThreads.append( newThread )
    newThread.start()
def set_keydown_handler( handler ):
    def inner_handler( e ):
        global no_current_keyhandler_call
        if no_current_keyhandler_call:
            no_current_keyhandler_call = False
            handler( e.keysym )
            no_current_keyhandler_call = True
    _getCanvas()._canvas._root().bind( "<Any-KeyPress>", inner_handler )
    _getCanvas()._canvas._root().update()
def unset_keydown_handler():
    _getCanvas()._canvas._root().unbind( "<Any-KeyPress>" )
def set_mousedown_handler( handler ):
    def inner_handler( e ):
        global no_current_mousehandler_call
        if no_current_mousehandler_call:
            no_current_mousehandler_call = False
            handler( e.x, e.y, e.num )
            no_current_mousehandler_call = True
    _getCanvas()._canvas.bind( "<Any-Button>", inner_handler )
    _getCanvas()._canvas._root().update()
def set_mouseup_handler( handler ):
    def inner_handler( e ):
        global no_current_mousehandler_call
        if no_current_mousehandler_call:
            no_current_mousehandler_call = False
            handler( e.x, e.y, e.num )
            no_current_mousehandler_call = True
    _getCanvas()._canvas.bind( "<Any-ButtonRelease>", inner_handler )
    _getCanvas()._canvas._root().update()
def set_mousemotion_handler( handler ):
    def inner_handler( e ):
        global no_current_mousehandler_call
        if no_current_mousehandler_call:
            no_current_mousehandler_call = False
            handler( e.x, e.y )
            no_current_mousehandler_call = True
    _getCanvas()._canvas.bind( "<Motion>", inner_handler )
    _getCanvas()._canvas._root().update()
#def (  ):      _getNavigator()._( )
