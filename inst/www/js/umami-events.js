/**
 * Umami Event Tracking Helper
 * Provides utility functions for tracking user interactions with Umami analytics
 */

// Wait for Umami to initialize before defining tracking functions
window.umamiReady = function() {
  return typeof window.umami !== 'undefined' && window.umami.track;
};

// Tracking function wrapper. The server calls these functions whether or not
// analytics are on (dashboard_analytics()), so they are always defined; without
// the Umami script on the page they do nothing (R22.18: they used to poll for
// it for 5 seconds, and log every event to the console).
window.trackEvent = function(eventName, properties) {
  if (window.umamiReady()) {
    try {
      window.umami.track(eventName, properties || {});
    } catch (e) {
      console.error('[Umami] Error tracking event:', eventName, e);
    }
  }
};

// Specialized tracking functions
//
// CODE_REVIEW.md R20.D12: no event carries a study string. Column names,
// treatment and reference-arm names and error messages (which can quote a
// column name) are unpublished study details; they are dropped here, the one
// place every event passes through, whatever the caller passes. Events carry
// only the analysis type, the method, counts, booleans and fixed codes.

/**
 * Track data file upload
 * @param {string} format - 'csv' or 'rda'
 * @param {number} rows - number of rows in dataset
 * @param {number} cols - number of columns
 */
window.trackDataUpload = function(format, rows, cols) {
  window.trackEvent('data_upload', {
    format: format,
    rows: rows,
    columns: cols,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track column mapping completion
 * @param {string} treatment - treatment column name
 * @param {string} time - time column name
 * @param {string} volume - volume column name
 * @param {boolean} survivalMapped - whether survival event column is mapped
 */
window.trackColumnMapping = function(treatment, time, volume, survivalMapped) {
  // R20.D12: which mapping options were used, never the column names.
  window.trackEvent('column_mapping_complete', {
    time_from_date: time === 'date',
    volume_calculated: volume === 'calculated',
    survival_event_mapped: survivalMapped === true,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track tab navigation
 * @param {string} tabName - name of tab being navigated to
 */
window.trackTabNavigation = function(tabName) {
  window.trackEvent('tab_navigate', {
    tab: tabName,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track analysis start
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {string} method - 'lme4', 'auc', 'cox', etc.
 * @param {string} referenceGroup - reference group name
 */
window.trackAnalysisStart = function(analysisType, method, referenceGroup) {
  window.analysisStartTime = Date.now(); // Store start time for duration calculation
  // R20.D12: the reference arm's name is never sent.
  window.trackEvent('analysis_run', {
    analysis_type: analysisType,
    method: method,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track analysis completion with runtime
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {boolean} success - whether analysis completed successfully
 * @param {string} errorMessage - error message if failed (optional)
 * @param {number} durationSeconds - optional override for duration (otherwise uses start time)
 */
window.trackAnalysisComplete = function(analysisType, success, errorMessage, durationSeconds) {
  var duration = durationSeconds;
  if (!duration && window.analysisStartTime) {
    duration = (Date.now() - window.analysisStartTime) / 1000;
    window.analysisStartTime = null;
  }
  
  var properties = {
    analysis_type: analysisType,
    success: success,
    duration_seconds: duration || 0,
    timestamp: new Date().toISOString()
  };
  
  // R20.D12: whether it failed, not the message.
  if (errorMessage) {
    properties.error = true;
  }
  
  window.trackEvent('analysis_complete', properties);
};

/**
 * Track download action
 * @param {string} analysisType - 'tumor_growth', 'survival', etc.
 * @param {string} contentType - 'plot', 'table', 'results', etc.
 * @param {string} format - 'png', 'pdf', 'csv', 'txt', etc.
 */
window.trackDownload = function(analysisType, contentType, format) {
  window.trackEvent('download_results', {
    analysis_type: analysisType,
    content_type: contentType,
    format: format,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track error condition
 * @param {string} errorType - 'validation', 'analysis', 'file_read', etc.
 * @param {string} errorMessage - description of error
 * @param {string} context - module/section where error occurred
 */
window.trackError = function(errorType, errorMessage, context) {
  // R20.D12: the message can quote a column name; only its type and the
  // module are sent.
  window.trackEvent('error_encountered', {
    error_type: errorType,
    context: context,
    timestamp: new Date().toISOString()
  });
};

/**
 * Track demo data loading
 * @param {string} demoName - name of demo dataset
 */
window.trackDemoDataLoaded = function(demoName) {
  window.trackEvent('demo_data_loaded', {
    demo_name: demoName,
    timestamp: new Date().toISOString()
  });
};
