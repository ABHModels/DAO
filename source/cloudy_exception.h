#ifndef CLOUDY_EXCEPTION_H
#define CLOUDY_EXCEPTION_H

// Catch-all macro for Cloudy exceptions.
// Usage:  try { ... cdEXIT(exit_status); } CLOUDY_CATCH_ALL(exit_status);
//
// This must appear after a try block that ends with cdEXIT().

#define CLOUDY_CATCH_ALL(exit_status) \
	catch( bad_alloc ) \
	{ \
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer " \
			"ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" ); \
		exit_status = ES_BAD_ALLOC; \
	} \
	catch( out_of_range& e ) \
	{ \
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n", \
			e.what() ); \
		exit_status = ES_OUT_OF_RANGE; \
	} \
	catch( domain_error& e ) \
	{ \
		fprintf( ioQQQ, " DISASTER - A vectorized math routine threw a domain_error. Bailing out...\n" ); \
		fprintf( ioQQQ, " What() = %s", e.what() ); \
		exit_status = ES_DOMAIN_ERROR; \
	} \
	catch( bad_assert& e ) \
	{ \
		MyAssert( e.file(), e.line(), e.comment() ); \
		exit_status = ES_BAD_ASSERT; \
	} \
	catch( bad_signal& e ) \
	{ \
		if( e.sig() == SIGILL ) \
		{ \
			if( ioQQQ != NULL ) \
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" ); \
			exit_status = ES_ILLEGAL_INSTRUCTION; \
		} \
		else if( e.sig() == SIGFPE ) \
		{ \
			if( ioQQQ != NULL ) \
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" ); \
			exit_status = ES_FP_EXCEPTION; \
		} \
		else if( e.sig() == SIGSEGV ) \
		{ \
			if( ioQQQ != NULL ) \
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" ); \
			exit_status = ES_SEGFAULT; \
		} \
		CLOUDY_CATCH_SIGBUS_(exit_status, e) \
		else \
		{ \
			if( ioQQQ != NULL ) \
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() ); \
			exit_status = ES_UNKNOWN_SIGNAL; \
		} \
	} \
	catch( cloudy_abort& e ) \
	{ \
		fprintf( ioQQQ, " ABORT DISASTER PROBLEM - Cloudy aborted, reason: %s\n", e.comment() ); \
		exit_status = ES_CLOUDY_ABORT; \
	} \
	catch( cloudy_exit& e ) \
	{ \
		if( ioQQQ != NULL ) \
		{ \
			ostringstream oss; \
			oss << " [Stop in " << e.routine(); \
			oss << " at " << e.file() << ":" << e.line(); \
			if( e.exit_status() == 0 ) \
				oss << ", Cloudy exited OK]"; \
			else \
				oss << ", something went wrong]"; \
			fprintf( ioQQQ, "%s\n", oss.str().c_str() ); \
		} \
		exit_status = e.exit_status(); \
	} \
	catch( std::exception& e ) \
	{ \
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n", \
			e.what() ); \
		exit_status = ES_UNKNOWN_EXCEPTION; \
	} \
	catch( ... ) \
	{ \
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" ); \
		exit_status = ES_UNKNOWN_EXCEPTION; \
	}

// Helper macro for SIGBUS (not available on all platforms)
#ifdef SIGBUS
#define CLOUDY_CATCH_SIGBUS_(exit_status, e) \
	else if( e.sig() == SIGBUS ) \
	{ \
		if( ioQQQ != NULL ) \
			fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" ); \
		exit_status = ES_BUS_ERROR; \
	}
#else
#define CLOUDY_CATCH_SIGBUS_(exit_status, e)
#endif

#endif // CLOUDY_EXCEPTION_H