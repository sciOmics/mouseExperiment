// Navbar layout and tab management helpers
// Extracted from inline JS in app.R (issue 3.3)

$(document).ready(function() {
  // Move Log, Help and About to a right-aligned nav list
  var $mainNav = $('.navbar-nav:first');
  var logLi      = $mainNav.find('a[data-value="Log"]').closest('li').detach();
  var helpLi     = $mainNav.find('a[data-value="Help"]').closest('li').detach();
  var citationLi = $mainNav.find('a[data-value="About"]').closest('li').detach();
  var $rightNav  = $('<ul class="nav navbar-nav navbar-right" id="nav-right"></ul>');
  $rightNav.append(logLi).append(helpLi).append(citationLi);
  // Bootstrap 5: no .navbar-collapse wrapper — append directly to .container-fluid
  var $container = $mainNav.closest('.container-fluid, .navbar-collapse');
  if ($container.length === 0) $container = $mainNav.parent();
  $container.append($rightNav);


  // Mark main nav for centering via CSS
  $mainNav.addClass('navbar-center');

  // Cross-nav active-state sync:
  // Bootstrap only clears 'active' within the same <ul>. Because Help/Citation
  // live in #nav-right (a separate <ul>), clicking a main-nav tab leaves them
  // visually active, causing a double-click requirement.
  // Solution: clear active from the *other* nav whenever a tab is clicked.
  $(document).on('click', '#nav-right a[data-value]', function() {
    $mainNav.find('li.active').removeClass('active');
  });
  $(document).on('click', '.navbar-center a[data-value]', function() {
    $('#nav-right li.active').removeClass('active');
  });
});

// Prevent clicks on disabled tabs
$(document).on('click', '.nav-tab-disabled > a, .nav-tab-disabled > .nav-link', function(e) {
  e.preventDefault(); e.stopPropagation(); return false;
});

function disableNavTab(val) {
  $('.navbar-nav li a[data-value="' + val + '"]').closest('li').addClass('nav-tab-disabled');
  $('#nav-right li a[data-value="' + val + '"]').closest('li').addClass('nav-tab-disabled');
}
function enableNavTab(val) {
  $('.navbar-nav li a[data-value="' + val + '"]').closest('li').removeClass('nav-tab-disabled');
  $('#nav-right li a[data-value="' + val + '"]').closest('li').removeClass('nav-tab-disabled');
}
window.disableNavTab = disableNavTab;
window.enableNavTab  = enableNavTab;
