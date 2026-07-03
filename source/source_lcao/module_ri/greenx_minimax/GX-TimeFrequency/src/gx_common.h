
#ifdef HAVE_FC_LONG_LINES
#define _REGISTER_EXC(msg) call register_exc(msg, filename=__FILE__, lineno=__LINE__) 
#else
#define _REGISTER_EXC(msg) call register_exc(msg)
#endif
