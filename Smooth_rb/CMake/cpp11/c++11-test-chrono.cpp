#include <chrono>

int main()
{
	uint64_t const start_time( std::time( nullptr ) );
	return ( start_time >= static_cast< uint64_t >( 0 ) ) ? 0 : 1;
}