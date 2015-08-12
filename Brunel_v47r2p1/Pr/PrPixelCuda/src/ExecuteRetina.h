/* 
 * File:   ExecuteRetina.h
 * Author: towelenee
 *
 * Created on March 15, 2015, 2:15 AM
 */
#pragma once
#include <stdint.h>
#include <vector>


int independent_execute(
    const std::vector<std::vector<uint8_t> > & input,
    std::vector<std::vector<uint8_t> > & output);

void independent_post_execute(const std::vector<std::vector<uint8_t> > & output);

int retina(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output);

/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input  
 * @param output 
 */
int retina_invocation(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output);
